% for testing the helper functions to do color to gray transformation

addpath(genpath('LE-Fusion'))

%% prep cifar 

imageDB = "cifar"; imNumber = 2;
% image 1 is a cat, image 2 is a bird
% imageDB = "local"; imNumber = "seahorse";
targetMode = "color";
sourceMode = "bw";
repType = "average";
repPoints = "TL";
% repPoints = ["TL", "TR", "BL", "BR"];
patchSize = [2 2];

spatialSigma = [1/5, 1/5];
% spatialSigma = [1/10, 1/10];
% spatialSigma = [0 0];
% useFusedImage = {true, "bicubic"};
useFusedImage = false;
maskFirst = false;
[X, Y, Z, images, idx] = prepImageData(imageDB, imNumber, targetMode, sourceMode,...
    patchSize, repType, repPoints, spatialSigma, useFusedImage, maskFirst);
% neighborSize = 100;
neighborSize = 40;
nodims = 15;
sigS = 0;
sigT = 0;
usedClust = false;
usedOOSFrameSubset = false;
%% look at images (bw - color)

% figure, imshow(images.rep, 'InitialMagnification', 'fit'), title('Low Res. Color')
% figure, imshow(images.source, 'InitialMagnification', 'fit'), title('High Res. BW Image')
% figure, imshow(images.target, 'InitialMagnification', 'fit'), title('High Res. Color')



%% create embedding
% sig = 0;
% nodims = 15;

mapTarget.X = X;
mapTarget.k = neighborSize;
mapTarget.sigma = 1.0*sigT;
% mapTarget.sigma = 0.14114;
% mapTarget.sigma = 0.026355;
mapTarget.name = 'Laplacian';
mapSource.X = Z;
mapSource.k = neighborSize;
mapSource.sigma = 1.0*sigS;
% mapSource.sigma = 0.083212;
% mapSource.sigma = 0.030351;
mapSource.name = 'Laplacian';

n = prod(size(images.target, [1 2]));

[mapSource.vec, mapSource.val, ~, mapSource.sigma, ~, mapSource.aff] = ...
    lapEig(mapSource.X, nodims, mapSource.k, 'sigma', mapSource.sigma, 'beQuiet', true);
[mapTarget.vec, mapTarget.val, ~, mapTarget.sigma, ~, mapTarget.aff] = ...
    lapEig(mapTarget.X, nodims, mapTarget.k, 'sigma', mapTarget.sigma, 'beQuiet', true);
mapSource.K = mapSource.aff;
mapSource.max_dist = sqrt(-(2*mapSource.sigma^2)*log(min(mapSource.K(mapSource.K ~= 0), [], 'all')));

%%
cN = 4;
if usedOOSFrameSubset
    mapSource.SampleSet = idx.maskSmall;
else
    mapSource.SampleSet = ones(size(X, 1), 1, "logical");
end
rotateMethod = "procrustes";
% rotateMethod = "scalerot";
clustMethod = "k-means";
% clustMethod = "curve";
% [mapSource, lapCustomS, idxS, idxT, cSInd, cTInd, proTG] = clustAndRotate(mapSource, mapTarget, Y, clustNum, clustType)
[mapSource, lapCustomS, idxS, idxT, cSInd, cTInd, proTG] = ...
    clustAndMatch(mapSource, mapTarget, Y, cN, clustMethod, rotateMethod);
usedClust = true;

%% 

preType = "matching";
usedClust = false;
%%
preType = "clustering";

%% preimages

n = prod(size(images.target, [1 2]));
pstart = tic;
PSy = zeros(n, size(X, 2));
% isfield(idx, "mask")
if preType == "clustering" && ~isfield(idx, "mask")
    PSy(idx.nz, :) = LEPre_init(mapSource.rotX, mapTarget);
    txtRot = clustMethod + " clustering";
elseif preType == "clustering" && isfield(idx, "mask")
    PSy(idx.nz & idx.mask, :) = LEPre_init(mapSource.rotX(idx.mask(idx.nz), :), mapTarget);
    txtRot = clustMethod + " clustering";
elseif preType == "matching" && isfield(idx, "mask") % default action is to only preimage mask points if available
    PSy(idx.nz & idx.mask, :) = LEPre_init(lapCustomS(idx.mask(idx.nz), :), mapTarget);
    txtRot = preType;
    clustMethod = "matching";
elseif preType == "matching" && ~isfield(idx, "mask")
    PSy(idx.nz, :) = LEPre_init(lapCustomS, mapTarget);
    txtRot = preType;
    clustMethod = "matching";
end
PSnz = PSy(idx.nz, :);
% PSy = LEPre_init(lapCustomg, mapC);
PSy = PSy(:,1:size(images.target, 3)); % cut off spatial dimensions
tPS = toc(pstart);
errorPS = norm(images.hX - PSy, 'fro')/norm(images.hX, 'fro');
usedPre = true;
%%
PSg = reshape(PSy, size(images.target));

if size(images.target, 3) == 3
    simPS = ssim(PSg, rgb2hsv(double(images.target)));
    psnrPS = psnr(PSg, rgb2hsv(double(images.target)));
else
    simPS = ssim(PSg, double(images.target));
    psnrPS = psnr(PSg, double(images.target));
end

badMask = any(PSg > 1, 3) | any(PSg < 0, 3) | any(isnan(PSg), 3);
badMask = ~badMask(:);
filterErrorPS = norm(images.hX(badMask, :) - PSy(badMask, :), 'fro')/norm(images.hX(badMask, :), 'fro');

% msePS = ((errorPS)*(norm(images.hX, 'fro')).^2 / numel(images.hX));
msePS = immse(double(images.hX), PSy);
% psnrPS = 20*log10(max(PSy, [], 'all')) -...
%     10*log10(msePS);


fprintf('Finished PS matching in %.3f sec. \nBad pixel count is %d. \n', tPS, n - nnz(badMask))
nx = size(images.source, 1);
ny = size(images.source, 2);
n = nx*ny;


%% display preimages
useSpatial = (size(images.hX, 2) < size(mapTarget.X, 2));

if ~useSpatial
    txtSpatial = "w/o spatial info";
    fileSpatial = "";
else
    txtSpatial = "with spatial info";
    fileSpatial = "Spatial";
end



sizevec = [reshape((repmat((1:ny)', [1, nx]))', [], 1), reshape(repmat((1:nx)', [1, ny]), [], 1)];


if isnan(errorPS) || isnan(simPS)
    txtError = "Error";
    txtSim = "Error";
    txtMse = "Error";
    txtPsnr = "Error";
else
    txtError = string(errorPS);
    txtSim = string(simPS);
    txtMse = string(msePS);
    txtPsnr = string(psnrPS);
end

badMask = any(PSg > 1, 3) | any(PSg < 0, 3) | any(isnan(PSg), 3);
badMask = ~badMask(:);
filterErrorPS = norm(images.hX(badMask, :) - PSy(badMask, :), 'fro')/norm(images.hX(badMask, :), 'fro');



isColor = size(images.hX, 2) > 1;
bigIndex = (1:n)';


dtPre = [dataTipTextRow("Pixel Location",literalNum(sizevec)), dataTipTextRow("Source Value",double(images.source(:))),...
         dataTipTextRow("True Target Value", literalNum(images.hX, 2)), dataTipTextRow("Computed Target Value", literalNum(PSy, 2)),...
         dataTipTextRow("Pixel Index", bigIndex), dataTipTextRow("Cluster #", idxS)];
if isColor
    dtPre(end + 1: end + 2) = [dataTipTextRow("True RGB", literalNum(hsv2rgb(images.hX), 3)), dataTipTextRow("Computed RGB", literalNum(hsv2rgb(PSg), 3))];
end


fG = figure('WindowState', 'maximized'); 
dcmG = datacursormode(fG);
if isColor
    zG = imshow(hsv2rgb(PSg), 'InitialMagnification', 'fit'); 
else
    zG = imshow(PSg, 'InitialMagnification', 'fit'); 
end
title(["Recovered hi-res " +  targetMode + " ( " + txtSpatial + " )",...
    "using " + txtRot + " rotation + shift",...
    "PSNR = " + txtPsnr + ", SSI = " + txtSim])
% title(["Recovered hi-res color (" + txtSpatial + ")",...
%     "using " + "matching"])
txtCaptionG = {'Relative Error = ' + txtError,'Error w/o Bad Pixels = ' + string(filterErrorPS), string(n - nnz(badMask)) + ' bad pixels',...
    'SS Index = ' + txtSim, "Mean Squared Error = " + txtMse,...
    "PSNR = " + txtPsnr + "dB", '{\it\sigma_G} = ' +  string(mapSource.sigma), '{\it\sigma_C} = ' +  string(mapTarget.sigma),...
     'NN = ' + string(mapSource.k), 'Dims = ' + string(nodims), 'NumClusters = '  + string(cN),...
      'PreImage Time = ' + string(tPS) + 's'};
% ,...
%       '{\it\mu_G} = ' +  string(muY*mapG.sigma) + ' (' + string(muY) + ')', '{\it\mu_C} = ' +  string(muX*mapC.sigma) + ...
%       ' (' + string(muY) + ')'};
% annotation(fG, 'textbox',[.05 .2 .2 .5], 'String', txtCaptionG, 'FitBoxToText', 'on')
dcmG.Enable = 'off'; datatip(zG); zG.DataTipTemplate.DataTipRows = dtPre; 
delete(findobj(gca(),'Type','DataTip'))


eG = figure; edG = datacursormode(eG);
if isColor
    ezG = imshow(rgb2gray(abs(hsv2rgb(double(images.target)) - hsv2rgb(PSg))), 'InitialMagnification', 'fit');
else
    ezG = imshow(abs(double(images.target) - PSg), 'InitialMagnification', 'fit');
end
title(["Error Image of " + sourceMode  + " Preimage (" + txtSpatial + ")", ...
    "using " + txtRot + " rotation + shift"])
% title(["Error Image of BW Preimage (" + txtSpatial + ")", ...
%     "using " + "matching"])
edG.Enable = 'off'; datatip(ezG); ezG.DataTipTemplate.DataTipRows = dtPre;
delete(findobj(gca(),'Type','DataTip'))


% %% color conv. filtering
% 
% fmat = ones(3); fmat(2, 2) = 0;
% patchCt = conv2(double(reshape(badMask, nx, ny)), rot90(fmat, 2), 'same');
% 
% filterPSg = convn(double(reshape(badMask, nx, ny)).*hsv2rgb(PSg), rot90(fmat, 2), 'same');
% filterPSg = filterPSg .* (1./patchCt);
% 
% filterPSg = rgb2hsv(filterPSg);
% chk = repmat(patchCt, [1 1 3]);
% chk = ~logical(chk);
% filterPSg(chk) = PSg(chk);
% filterPSg = double(reshape(badMask, nx, ny)).*PSg +...
%     double(reshape(~badMask, nx, ny)).*filterPSg;
% 
% 
% disp(psnrPS)
% disp(simPS)
% if size(images.target, 3) == 3
%     filterSimPS = ssim(filterPSg, rgb2hsv(double(images.target)))
%     filterPsnrPS = psnr(filterPSg, rgb2hsv(double(images.target)))
% else
%     filterSimPS = ssim(filterPSg, double(images.target))
%     filterPsnrPS = psnr(filterPSg, double(images.target))
% end
% 
% ffil = figure; imshow(hsv2rgb(filterPSg), 'InitialMagnification', 'fit')
% title(["Recovered hi-res " +  targetMode + " ( " + txtSpatial + " )",...
%     "using " + txtRot + " rotation + shift " + "+ 3x3 average filtering",...
%     "PSNR = " + string(filterPsnrPS) + ", SSI = " + string(filterSimPS)])
% %% replace bad pixel with value from low res representative
% 
% % simple resizing
% replacePSg = zeros(size(images.target));
% % A = zeros(size(images.target));
% for i = 1:size(images.rep, 3)
%     replacePSg(:,:,i) = kron(images.rep(:,:,i), ones(2));
% end
% 
% replacePSg = rgb2hsv(replacePSg);
% replacePSg = double(reshape(badMask, nx, ny)).*PSg +...
%     double(reshape(~badMask, nx, ny)).*replacePSg;
% 
% 
% disp(psnrPS)
% disp(simPS)
% if size(images.target, 3) == 3
%     replaceSimPS = ssim(replacePSg, rgb2hsv(double(images.target)))
%     replacePsnrPS = psnr(replacePSg, rgb2hsv(double(images.target)))
% else
%     replaceSimPS = ssim(replacePSg, double(images.target))
%     replacePsnrPS = psnr(replacePSg, double(images.target))
% end
% 
% figure, imshow(hsv2rgb(replacePSg), 'InitialMagnification', 'fit')
% title(["Recovered hi-res " +  targetMode + " ( " + txtSpatial + " )",...
%     "using " + txtRot + " rotation + shift ", "Bad pixels replaced with Low Res. Value",...
%     "PSNR = " + string(replacePsnrPS) + ", SSI = " + string(replaceSimPS)])
% 
% %% combine both approaches, replace bad pixel with average including low res value
% fmat = ones(3); fmat(2, 2) = 0;
% lowPatchCt = conv2(double(reshape(badMask, nx, ny)), rot90(fmat, 2), 'same');
% lowPatchCt = lowPatchCt + 1;
% avgLowPSg = zeros(size(images.target));
% for i = 1:size(images.rep, 3)
%     avgLowPSg(:,:,i) = kron(images.rep(:,:,i), ones(2));
% end
% 
% avgLowPSg = avgLowPSg + convn(double(reshape(badMask, nx, ny)).*hsv2rgb(PSg), rot90(fmat, 2), 'same');
% avgLowPSg = avgLowPSg .* (1./lowPatchCt);
% 
% avgLowPSg = rgb2hsv(avgLowPSg);
% avgLowPSg = double(reshape(badMask, nx, ny)).*PSg +...
%     double(reshape(~badMask, nx, ny)).*avgLowPSg;
% 
% 
% disp(psnrPS)
% disp(simPS)
% if size(images.target, 3) == 3
%     avgLowSimPS = ssim(avgLowPSg, rgb2hsv(double(images.target)))
%     avgLowPsnrPS = psnr(avgLowPSg, rgb2hsv(double(images.target)))
% else
%     avgLowSimPS = ssim(avgLowPSg, double(images.target))
%     avgLowPsnrPS = psnr(avgLowPSg, double(images.target))
% end
% 
% figure, imshow(hsv2rgb(replacePSg), 'InitialMagnification', 'fit')
% title(["Recovered hi-res " +  targetMode + " ( " + txtSpatial + " )",...
%     "using " + txtRot + " rotation + shift ", "Bad pixels replaced with 3x3 avg + Low Res. Value",...
%     "PSNR = " + string(avgLowPsnrPS) + ", SSI = " + string(avgLowSimPS)])


%% set up plot variables 

% if (usedClust && usedClustColor)
%     warning('Make sure the correct clustering is used. \n')
% end
imDim = size(images.rep, [1 2]);
sizevec = [reshape((repmat((1:ny)', [1, nx]))', [], 1), reshape(repmat((1:nx)', [1, ny]), [], 1)];
sizevecColor = [reshape((repmat((1:imDim(2))', [1, imDim(1)]))', [], 1), reshape(repmat((1:imDim(1))', [1, imDim(2)]), [], 1)];

if ~(usedClust)
    cN = 1; % default to 1 cluster if no clustering is used
    idxS = true(size(Y, 1), 1); cSInd = idxS;
    idxT = true(size(Z, 1), 1);  cTInd = idxT;
    cNames = {'c1'};
    txtRot = "not clustered ";
end

clustView = 1:cN; % which clusters to plot?

if ~usedPre
    PSy = zeros(size(images.hX));
end

if ~useSpatial
    txtSpatial = "w/o spatial info";
    fileSpatial = "Spectral";
else
    txtSpatial = "with spatial info";
    fileSpatial = "Spatial";
end


% coloring the different sized sets
redColors = ones(nnz(idx.nz), 3);
blueColors = ones(size(mapTarget.X, 1), 3);
for nn = 1:cN
    redColors(cSInd(:, nn), 1) = ((nn - 1)/(6*cN))*ones(nnz(cSInd(:, nn)), 1);
    blueColors(cTInd(:, nn), 1) = (2/3) - ((nn - 1)/(6*cN))*ones(nnz(cTInd(:, nn)), 1);
end
% redColors(bmidx, :) = zeros(numel(midx), 3); % rep points in black
% rep points in black separately
redColors = hsv2rgb(redColors);
blueColors = hsv2rgb(blueColors);

vD = [1, 2, 3];

%% look at graphs



p1 = struct(sourceMode, [], targetMode, []);
f1 = figure('WindowState','maximized'); 
for jj = 1:numel(clustView)
    nn = clustView(jj);
    gi = cSInd(:, nn); % temp for current source index (in largest set)
    ci = cTInd(:, nn); % temp for current target index (in largest set)
     
    

    dtRot = [dataTipTextRow("Pixel X-Location",sizevec(gi,1)),dataTipTextRow("Pixel Y-Location",sizevec(gi,2)),...
        dataTipTextRow(sourceMode + " Index",Y(gi, :)), dataTipTextRow("Cluster #", idxS(gi)),...
        dataTipTextRow("Computed " + targetMode + " Index", PSnz(gi,1))];
    dtRotColor = [dataTipTextRow("Pixel X-Location",sizevecColor(ci,1)),dataTipTextRow("Pixel Y-Location", sizevecColor(ci,2)),...
        dataTipTextRow("Cluster #", idxT(ci)),...
        ];



    p1.(sourceMode).(mapSource.cNames{nn}) = scatter3(mapSource.vecX(gi, vD(1)), mapSource.vecX(gi, vD(2)), mapSource.vecX(gi, vD(3)),...
        6, redColors(gi, :), 'o', 'DisplayName', mapSource.cNames{nn});
    hold on

    p1.(targetMode).(mapSource.cNames{nn}) = scatter3(mapTarget.vec(ci, vD(1)), mapTarget.vec(ci, vD(2)), mapTarget.vec(ci, vD(3)),...
        6, blueColors(ci, :), 'o', 'DisplayName', "Color " + mapSource.cNames{nn});
    hold on
    p1.rep.(mapSource.cNames{nn}) = scatter3(mapSource.vec(ci, vD(1)), mapSource.vec(ci, vD(2)), mapSource.vec(ci, vD(3)),...
        6, 'k', 'o');
    p1.(sourceMode).(mapSource.cNames{nn}).DataTipTemplate.DataTipRows = dtRot;
    p1.(targetMode).(mapSource.cNames{nn}).DataTipTemplate.DataTipRows = dtRotColor;
    p1.rep.(mapSource.cNames{nn}).DataTipTemplate.DataTipRows = dtRotColor;
   
end
axis equal
% legend show
title([sourceMode + " pixels"...
    "(" + string(size(mapSource.vecX, 1)) + ", rep. ( " + repType +  " ) pixels black, others red, points in red out-of-sample)",...
    " " + targetMode + " pixels (" + string(size(Z, 1)) + ", blue)", ...
    txtSpatial])
hold off

% % for printing:
% svPath = 'C:/Users/Jerry/Desktop/feb2Examples/feb2_';
% figure('visible', 'off'), 
% for jj = 1:numel(clustView)
%     nn = clustView(jj);
%     gi = cGInd(:, nn); % temp for current gray index 
%     ci = cCInd(:, nn); % temp for current color index 
% 
%     dtRot = [dataTipTextRow("Pixel X-Location",sizevec(gi,1)),dataTipTextRow("Pixel Y-Location",sizevec(gi,2)),...
%         dataTipTextRow("Gray Value",double(Ig(gi))), dataTipTextRow("Cluster #", idxG(gi)),...
%         dataTipTextRow("Computed Hue", PSy(gi,1)), dataTipTextRow("Computed Saturation",PSy(gi,2)),...
%         dataTipTextRow("Computed Value",PSy(gi,3))];
%     dtRotColor = [dataTipTextRow("Pixel X-Location",sizevecColor(ci,1)),dataTipTextRow("Pixel Y-Location", sizevecColor(ci,2)),...
%         dataTipTextRow("Cluster #", idxC(ci)),...
%         dataTipTextRow("Hue", X(ci,1)), dataTipTextRow("Saturation",X(ci,2)), dataTipTextRow("Value",X(ci,3)), ...
%         dataTipTextRow("Computed Hue", PSy(gi & bmidx,1)), dataTipTextRow("Computed Saturation",PSy(gi & bmidx,2)),...
%         dataTipTextRow("Computed Value",PSy(gi & bmidx,3))];
% 
%     p1.gray.(cNames{nn}) = scatter3(mapG.vecX(gi, vD(1)), mapG.vecX(gi, vD(2)), mapG.vecX(gi, vD(3)),...
%         6, redColors(gi, :), 'o', 'DisplayName', cNames{nn});
%     hold on
% 
%     p1.color.(cNames{nn}) = scatter3(mapC.vec(ci, vD(1)), mapC.vec(ci, vD(2)), mapC.vec(ci, vD(3)),...
%         6, blueColors(ci, :), 'o', 'DisplayName', "Color " + cNames{nn});
%     p1.gray.(cNames{nn}).DataTipTemplate.DataTipRows = dtRot;
%     p1.color.(cNames{nn}).DataTipTemplate.DataTipRows = dtRotColor;
%    
% end
% axis equal
% % legend show
% title(["Gray pixels"...
%     "("+ string(numel(Ig)) + ", top-left pixels black, others red, points in red out-of-sample)",...
%     "color pixels (" + string(numel(Ic)/3) + ", blue)"])
% hold off
% print(gcf, '-djpeg', svPath + imName + '_BW_and_Color.jpg', '-r0')
% close all

% 
% 
if (usedClust)
    p2 = struct(sourceMode, [], targetMode, []);
    f2 = figure('WindowState','maximized'); 
    for jj = 1:numel(clustView)
        nn = clustView(jj);
        gi = cSInd(:, nn); % temp for current T2 index (in largest set)
        ci = cTInd(:, nn); % temp for current T1 index (in largest set)
         
        
    
        dtRot = [dataTipTextRow("Pixel X-Location",sizevec(gi,1)),dataTipTextRow("Pixel Y-Location",sizevec(gi,2)),...
            dataTipTextRow(sourceMode + " Index",Y(gi, :)), dataTipTextRow("Cluster #", idxS(gi)),...
            dataTipTextRow("Computed " + targetMode + " Index", PSnz(gi,1))];
        dtRotColor = [dataTipTextRow("Pixel X-Location",sizevecColor(ci,1)),dataTipTextRow("Pixel Y-Location", sizevecColor(ci,2)),...
            dataTipTextRow("Cluster #", idxT(ci)),...
            ];

        repRot = ( mapSource.vec(cTInd(:, nn), :) )*proTG.(mapSource.cNames{nn}).T +...
            proTG.(mapSource.cNames{nn}).c(1, :);
    
        p2.(sourceMode).(mapSource.cNames{nn}) = scatter3(mapSource.rotX(gi, vD(1)), mapSource.rotX(gi, vD(2)), mapSource.rotX(gi, vD(3)),...
            6, redColors(gi, :), 'o', 'DisplayName', mapSource.cNames{nn});
        hold on
    
        p2.(targetMode).(mapSource.cNames{nn}) = scatter3(mapTarget.vec(ci, vD(1)), mapTarget.vec(ci, vD(2)), mapTarget.vec(ci, vD(3)),...
            6, blueColors(ci, :), 'o', 'DisplayName', "Color " + mapSource.cNames{nn});
        p2.rep.(mapSource.cNames{nn}) = scatter3(repRot(:, vD(1)), repRot(:, vD(2)), repRot(:, vD(3)),...
            6, 'k', 'o');
        p2.(sourceMode).(mapSource.cNames{nn}).DataTipTemplate.DataTipRows = dtRot;
        p2.(targetMode).(mapSource.cNames{nn}).DataTipTemplate.DataTipRows = dtRotColor;
        p2.rep.(mapSource.cNames{nn}).DataTipTemplate.DataTipRows = dtRotColor;
       
    end
    axis equal
%     legend show
    title([txtRot + " Rotated + Shifted " + sourceMode + " pixels"...
        "(" + string(size(mapSource.vecX, 1)) + ", rep. ( " + repType +  " ) pixels black, points in red out-of-sample)", ...
        " " + targetMode + " pixels (" + string(size(Z, 1)) + ", blue)", ...
        txtSpatial])
    hold off
    
%     % for printing:
%     svPath = 'C:/Users/Jerry/Desktop/feb2Examples/feb2_';
%     figure('visible', 'off'),
%     for jj = 1:numel(clustView)
%         nn = clustView(jj);
%         gi = cGInd(:, nn); % temp for current gray index 
%         ci = cCInd(:, nn); % temp for current color index 
%     
%         dtRot = [dataTipTextRow("Pixel X-Location",sizevec(gi,1)),dataTipTextRow("Pixel Y-Location",sizevec(gi,2)),...
%             dataTipTextRow("Gray Value",double(Ig(gi))), dataTipTextRow("Cluster #", idxG(gi)),...
%             dataTipTextRow("Computed Hue", PSy(gi,1)), dataTipTextRow("Computed Saturation",PSy(gi,2)),...
%             dataTipTextRow("Computed Value",PSy(gi,3))];
%         dtRotColor = [dataTipTextRow("Pixel X-Location",sizevecColor(ci,1)),dataTipTextRow("Pixel Y-Location", sizevecColor(ci,2)),...
%             dataTipTextRow("Cluster #", idxC(ci)),...
%             dataTipTextRow("Hue", X(ci,1)), dataTipTextRow("Saturation",X(ci,2)), dataTipTextRow("Value",X(ci,3)), ...
%             dataTipTextRow("Computed Hue", PSy(gi & bmidx,1)), dataTipTextRow("Computed Saturation",PSy(gi & bmidx,2)),...
%             dataTipTextRow("Computed Value",PSy(gi & bmidx,3))];
%     
%         p2.gray.(cNames{nn}) = scatter3(mapG.rotX(gi, vD(1)), mapG.rotX(gi, vD(2)), mapG.rotX(gi, vD(3)),...
%             6, redColors(gi, :), 'o', 'DisplayName', cNames{nn});
%         hold on
%     
%         p2.color.(cNames{nn}) = scatter3(mapC.vec(ci, vD(1)), mapC.vec(ci, vD(2)), mapC.vec(ci, vD(3)),...
%             6, blueColors(ci, :), 'o', 'DisplayName', "Color " + cNames{nn});
%         p2.gray.(cNames{nn}).DataTipTemplate.DataTipRows = dtRot;
%         p2.color.(cNames{nn}).DataTipTemplate.DataTipRows = dtRotColor;
%        
%     end
%     axis equal
% %     legend show
%     title(["MATLAB Rotated + Shifted Gray pixels"...
%         "("+ string(numel(Ig)) + ", top-left pixels in black, points in red out-of-sample)", ...
%         " color pixels (" + string(numel(Ic)/3) + ", blue)"])
%     hold off
%     print(gcf, '-djpeg', svPath + imName + '_BW_Rotated.jpg', '-r0')
%     close all
end
% 
% p3 = struct('gray', [], 'color', []);
p3 = struct(sourceMode, [], targetMode, []);
f3 = figure('WindowState','maximized'); 
for jj = 1:numel(clustView)
    nn = clustView(jj);
    gi = cSInd(:, nn); % temp for current T2 index (in largest set)
    ci = cTInd(:, nn); % temp for current T1 index (in largest set)
     
    

        dtRot = [dataTipTextRow("Pixel X-Location",sizevec(gi,1)),dataTipTextRow("Pixel Y-Location",sizevec(gi,2)),...
            dataTipTextRow(sourceMode + " Index",Y(gi, :)), dataTipTextRow("Cluster #", idxS(gi)),...
            dataTipTextRow("Computed " + targetMode + " Index", PSnz(gi,1))];
        dtRotColor = [dataTipTextRow("Pixel X-Location",sizevecColor(ci,1)),dataTipTextRow("Pixel Y-Location", sizevecColor(ci,2)),...
            dataTipTextRow("Cluster #", idxT(ci)),...
            ];

    p3.(sourceMode).(mapSource.cNames{nn}) = scatter3(lapCustomS(gi, vD(1)), lapCustomS(gi, vD(2)), lapCustomS(gi, vD(3)),...
        6, redColors(gi, :), 'o', 'DisplayName', mapSource.cNames{nn});
    hold on

    p3.(targetMode).(mapSource.cNames{nn}) = scatter3(mapTarget.vec(ci, vD(1)), mapTarget.vec(ci, vD(2)), mapTarget.vec(ci, vD(3)),...
        6, blueColors(ci, :), 'o', 'DisplayName', "Color " + mapSource.cNames{nn});
    p3.(sourceMode).(mapSource.cNames{nn}).DataTipTemplate.DataTipRows = dtRot;
    p3.(targetMode).(mapSource.cNames{nn}).DataTipTemplate.DataTipRows = dtRotColor;
   
end
axis equal
% legend show
title(["Perfect Matched " + sourceMode + " pixels"...
    "(" + string(size(mapSource.vecX, 1)) + ", rep. ( " + repType +  " ) pixels black, points in red out-of-sample)", ...
        " " + targetMode + " pixels (" + string(size(Z, 1)) + ", blue)", ...
        txtSpatial])
hold off

