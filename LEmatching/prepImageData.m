function [targetSet, sourceSet, repSet, images, index] = prepImageData(imageDB, imageName, targetMode, sourceMode, patchSize, repType, repPoints, spatialSigma, fusedImage, mskFirst)
% working on a suite of functions to do the tasks related to the image
% modality transformation problem

% need a function that provides an image and the sets X, Y, and Z to be
% processed. X is the lowres modality, Y is the high res modality, and Z is
% the set of representative points of Y (in correspondence to X). in
% addition, it needs to have the option of using spatial info

% images need to get output

% 
% inputs:
% global vars:
% patchSize, repType, repPoints, useSpatial, muX, muY,
% targetMode, sourceMode
% patchSize = [2, 2];
% neighborSize = 40;
% repType = "sample";
% repPoints = "TL";

% if isstring(repPoints)
%     repPoints = char(repPoints);
% end

if iscell(fusedImage)
    useFusedImage = fusedImage{1};
    interType = fusedImage{2};
else
    useFusedImage = fusedImage;
    interType = "bilinear";
end

useSpatial = true;
% useFusedImage = true;
muX = spatialSigma(1); muY = spatialSigma(2);
dbList = ["cifar", "local", "brats", "brainweb";...
        "color", "color", "mri", "mri"];
imageType = dictionary(dbList(1, :), dbList(2, :));
switch imageType(imageDB)
    case "color"
        [targetSet, sourceSet, repSet, images, index] = getColorImageData(imageDB, imageName);
    case "mri"
        if numel(imageName) < 2
            % use default volume = 1, slice = 68
            volNum = 1; sl = []; viewAx = 1;
        elseif numel(imageName) < 3
            volNum = imageName(1);
            sl = imageName(2); viewAx = 1;
        else
            volNum = imageName(1);
            sl = imageName(2); viewAx = imageName(3);
        end
        [targetSet, sourceSet, repSet, images, index] = getMRIData(imageDB, volNum, sl, viewAx, targetMode, sourceMode, repType, repPoints);
end
% if imageDB == "brats"
%     if numel(imageName) < 2
%         % use default volume = 1, slice = 68
%         volNum = 1; sl = 68;
%     else
%         volNum = imageName(1);
%         sl = imageName(2);
%     end
% 
%     [targetSet, sourceSet, repSet, images, index] = getMRIData(imageDB, volNum, sl, targetMode, sourceMode, repType, repPoints);
% else
%     [targetSet, sourceSet, repSet, images, index] = getColorImageData(imageDB, imageName);
% end

    function [target, source, rep, images, index] = getColorImageData(dset, imName)
    switch dset
        case "cifar"
            B = matfile('cifar_batch1.mat');
            kk = 1000;
            numImages = 100;
            rng(4+20+2022)
            im = randperm(10000, kk);
            iSet = randperm(kk, numImages); % a bit silly 
            Ic = B.dataBatch1(im(iSet(imName)), :); 
            Ic = permute(reshape(Ic, 32, 32, 3), [2 1 3]);
            [nx, ny] = size(Ic, 1:2);
            isSquare = @(x) size(x, 1) == size(x, 2);
            if isSquare(Ic)
                mm = nx;
            end
            imDim = ceil([nx, ny] ./ patchSize);
            if isSquare(Ic) && (patchSize(1) == patchSize(2))
                mc = imDim(1);
            end
            n = nx*ny;
            imName = "image" + string(imName);
        case "local"
            Ic = imread('horse_64.jpg');
            imName = "seahorse";
            % for subsampling
            patchSize = [2, 2];
            [nx, ny] = size(Ic, 1:2);
            % isSquare = nx == ny;
            isSquare = @(x) size(x, 1) == size(x, 2);
            if isSquare(Ic)
                mm = nx;
            end
            imDim = ceil([nx, ny] ./ patchSize);
            if isSquare(Ic) && (patchSize(1) == patchSize(2))
                mc = imDim(1);
            end
            n = nx*ny;
    end

    images = struct();
    images.imName = imName;
    if isinteger(Ic)
        hIc = double(Ic)./(pow2(8*(whos('Ic').bytes)/numel(Ic)) - 1);
    end
    hIc = rgb2hsv(double(hIc));
    images.hX = reshape(hIc, n, []);
    Ig = double(reshape(Ic, [], size(Ic, 3)));
    if size(Ic, 3) == 3
        Ig = reshape(Ig*[0.2989; 0.5870; 0.1140], size(Ic, 1), size(Ic, 2)); % need to start with a double to do the multiplication
    end
    Ig = Ig./255;
    images.source = Ig;
    images.target = hsv2rgb(hIc);
    bigIndex = (1:n)'; % for any miscellaneous indexing needs
    
    % X is points for target embedding, Y is full set of out-of-sample points,
    % Z is set of representative points
    index = struct();
    if repType == "sample"
        [midx, grpLabel] = patchLabel([nx, ny], patchSize, repPoints);
        if numel(midx) == (imDim(1)-1)*imDim(2)
            imDim = imDim - [1 0];
        elseif numel(midx) == imDim(1)*(imDim(2)-1)
            imDim = imDim - [0 1];
        elseif numel(midx) == (imDim(1) - 1)*(imDim(2)-1)
            imDim = imDim - 1;
        end
        Y = double(reshape(Ig, [], size(Ig, 3)));
        Ic = reshape(hIc, [], 3);
        Ic = Ic(midx, :);
        Ic = reshape(Ic, imDim(1), imDim(2), 3);
        images.rep = hsv2rgb(Ic);
        hX = double(reshape(hIc, [], 3));
        X = double(reshape(Ic, [], 3));
        if useFusedImage
            fY = imresize(Ic, [nx ny], interType, "Antialiasing", false);
            fY = reshape(fY, size(Y, 1), []);
            Y = double([Y, fY]);
        end
        bmidx = ismember(bigIndex, midx);
        % spatial info
        sizevec = [reshape((repmat((1:ny)', [1, nx]))', [], 1), reshape(repmat((1:nx)', [1, ny]), [], 1)];
        sizevecSmall = [reshape((repmat((1:imDim(2))', [1, imDim(1)]))', [], 1), reshape(repmat((1:imDim(1))', [1, imDim(2)]), [], 1)];
        
        muX = useSpatial*muX; muY = useSpatial*muY;
        % to turn off spatial info
        if muX == 0
            addX = zeros(prod(imDim), 0);
        else
            addX = muX.*sizevecSmall*diag(1./[imDim(2), imDim(1)]);
        end
        if muY == 0
            addY = zeros(n, 0);
        else
            addY = muY.*sizevec*diag(1./[ny, nx]);
        end
        target = [X, addX];
        source = [Y, addY];
        rep = [Y(midx, :), addY(midx, :)];
        index.rep = bmidx;
        index.nz = true(n, 1);
    elseif repType == "average"
        Y = double(reshape(Ig, [], size(Ig, 3)));
%         Ic = reshape(hIc, [], 3);
        Ic = imresize(hsv2rgb(hIc), imDim, 'bilinear', 'Antialiasing', false);
        images.rep = Ic;
        Ic = rgb2hsv(Ic);
        if useFusedImage
            fY = imresize(Ic, [nx ny], interType, "Antialiasing", false);
            rfY = reshape(fY, size(Y, 1), []);
            Y = double([Y, rfY]);
        end
        X = reshape(Ic, [], size(Ic, 3));
        if useFusedImage
            Z = imresize(cat(3, Ig, fY), imDim, 'bilinear', 'Antialiasing', false);
        else
            Z = imresize(Ig, imDim, 'bilinear', 'Antialiasing', false);
        end
        Z = reshape(Z, prod(imDim), []);
        sizevec = [reshape((repmat((1:ny)', [1, nx]))', [], 1), reshape(repmat((1:nx)', [1, ny]), [], 1)];
        sizevecSmall = [reshape((repmat((1:imDim(2))', [1, imDim(1)]))', [], 1), reshape(repmat((1:imDim(1))', [1, imDim(2)]), [], 1)];
        
        muX = useSpatial*muX; muY = useSpatial*muY;
        % to turn off spatial info
        % to turn off spatial info
        if muX == 0
            addX = zeros(prod(imDim), 0);
        else
            addX = reshape(sizevec, nx, ny, []);
            addX = imresize(addX, imDim, 'bilinear', 'Antialiasing', false);
            addX = muX.*reshape(addX, prod(imDim), [])*diag(1./[ny nx]');
        end
        if muY == 0
            addY = zeros(n, 0);
        else
            addY = muY.*sizevec./mm;
        end
        target = [X, addX];
        source = [Y, addY];
        rep = [Z, addX];
        index.rep = [];
        index.nz = true(n, 1);
    end
    

end
%%
% MRI

% inputs:
    function [target, source, rep, images, index] = getMRIData(imageDB, volumeNumber, slice, viewPlane, targetMode, sourceMode, repType, repPoints)
%     volumeNumber = 1;
%     slice = 68;
%     sourceMode = "T2w";
%     targetMode = "T1w";
%     repPoints = "TL";
    [~, ret] = system('hostname'); % for identification purposes
    ret = string(ret(1:end-1));
    
    if ret == "Jerky2"
        bratsDir = "C:/Users/Jerry/Documents/UMD Research/"; % to avoid cd problems
    elseif contains(ret, "theta")
        bratsDir = "/home/jemidih/MATLAB/";
    end
    switch imageDB
        case "brats"
            volPath = bratsDir + "Task01_BrainTumour/imagesTr/BRATS_" +...
                sprintf('%03d', volumeNumber) + ".nii.gz";
            I = niftiread(volPath);     
            if mskFirst
                kPath = bratsDir + "Task01_BrainTumour/CSFMask/BRATS_" +...
                    sprintf('%03d', volumeNumber) + "Mask.nii.gz";
                A = niftiread(kPath);
                I = I.*single(A);
            end
            I = permute(I, [circshift([1 2 3], 1 - viewPlane), 4]);
            planeShift = viewPlane < 3;
            I = permute(I, [circshift([1 2], -1*planeShift), 3 4]);
            I = flip(I, 1);
%             gtI = niftiread("C:\Users\Jerry\Documents\UMD Research\Task01_BrainTumour\labelsTr\BRATS_" + sprintf('%03d', volumeNumber) + ".nii.gz");
            
            if isempty(slice) || slice == 0
                slice = ceil(size(I, 3)/2);
            end
            I = I(:,:,slice,:);
            % gtI = permute(gtI(:,:,slice,:), [2 1 3 4]);
            for i = 1:size(I, 3)*size(I, 4)
                j = ceil(i/size(I, 4));
                jj = rem(i - 1, size(I, 4)) + 1;
                I(:,:,j,jj) = double(mat2gray(I(:,:,j,jj)));
            end
            mrModes = ["FLAIR", "T1w", "T1Gd", "T2w"];
            % the order is FLAIR, T1w, T1Gd, T2w via dataset.json file
            It1 = I(:,:,:,mrModes == targetMode);
            It2 = I(:,:,:,mrModes == sourceMode);
            imName = "slice" + string(slice);
            msk = niftiread(bratsDir + "Task01_BrainTumour/CSFMask/BRATS_" +...
                sprintf('%03d', volumeNumber) + "Mask.nii.gz");
            msk = permute(msk, [circshift([1 2 3], 1 - viewPlane), 4]);
            msk = permute(msk, [circshift([1 2], -1*planeShift), 3 4]);
            msk = flip(msk, 1);
            msk = msk(:,:,slice,1);
        case "brainweb"
            modePath = dictionary(["T1w", "T2w"],...
                ["brainweb/t1_icbm_normal_1mm_pn0_rf20.nii", "brainweb/t2_icbm_normal_1mm_pn0_rf20.nii"]);
            if volumeNumber >=4 % for MS volume
                modePath = dictionary(["T1w", "T2w"],...
                    ["brainweb/t1_ai_msles2_1mm_pn0_rf20.nii", "brainweb/t2_ai_msles2_1mm_pn0_rf20.nii"]);
            end
            mskPath = ["brainweb/icbm_normal_1mm_pn0_rf20BrainMaskver3Old.nii",...
                "brainweb/icbm_normal_1mm_pn0_rf20ThreshMask.nii",...
                "brainweb/icbm_normal_1mm_pn0_rf20BrainMaskver3.nii"];
            % skull strip mask at "brainweb/icbm_normal_1mm_pn0_rf20BrainMaskver3.nii"
%             fidT = fopen(modePath(targetMode));
%             It1 = fread(fidT, inf, 'uint8=>uint8');
% %             It1 = permute(reshape(It1(end:-1:1), [181, 217, 181]), [1 2 3]);
%             It1 = permute(reshape(It1(:), [181, 217, 181]), [1 2 3]);
            It1 = niftiread(modePath(targetMode));
            if volumeNumber >= 2
                msk = uint8(niftiread(mskPath(1)));
                It1 = msk.*It1;
            end
            % if volumeNumber == 3 % 
            if mskFirst
                It1 = niftiread(modePath(targetMode));
                msk = uint8(niftiread(mskPath(3)));
                It1 = msk.*It1;
                msk = uint8(niftiread(mskPath(2)));
                It1 = msk.*It1;
            end
            It1 = shiftdim(It1, viewPlane - 1);
            planeShift = viewPlane < 3;
            It1 = permute(It1, [circshift([1 2], -1*planeShift), 3]);
            It1 = flip(It1, 1);
%             It1 = flip(It1, 2);
            if isempty(slice) || slice == 0
                slice = ceil(size(It1, 3)/2);
            end
            It1 = It1(:, :, slice);
            It1(It1 < 0.1*max(It1(:))) = 0;
            if isinteger(It1)
                It1 = double(It1)./(pow2(8*(whos('It1').bytes)/numel(It1)) - 1);
            end
%             It1 = maskon(It1);
%             It1 = imresize(It1, 0.5, "bilinear");

%             fidS = fopen(modePath(sourceMode));
%             It2 = fread(fidS, inf, 'uint8=>uint8');
% %             It2 = permute(reshape(It2(end:-1:1), [181, 217, 181]), [1 2 3]);
%             It2 = permute(reshape(It2(:), [181, 217, 181]), [1 2 3]);
            It2 = niftiread(modePath(sourceMode));
            if volumeNumber >= 2
                msk = uint8(niftiread(mskPath(1)));
                It2 = msk.*It2;
            end
            % if volumeNumber == 3
            if mskFirst
                It2 = niftiread(modePath(sourceMode));
                msk = uint8(niftiread(mskPath(3)));
                It2 = msk.*It2;
                msk = uint8(niftiread(mskPath(2)));
                It2 = msk.*It2;
            end
            msk = uint8(niftiread(mskPath(2)));
            It2 = shiftdim(It2, viewPlane - 1);
            It2 = permute(It2, [circshift([1 2], -1*planeShift), 3]);
            It2 = flip(It2, 1);
%             It2 = flip(It2, 2);
            It2 = It2(:, :, slice);
            It2(It2 < 0.1*max(It2(:))) = 0;
            if isinteger(It2)
                It2 = double(It2)./(pow2(8*(whos('It2').bytes)/numel(It2)) - 1);
            end
            msk = uint8(niftiread(mskPath(3))) .* uint8(niftiread(mskPath(2)));
            msk = shiftdim(msk, viewPlane - 1);
            msk = permute(msk, [circshift([1 2], -1*planeShift), 3]);
            msk = flip(msk, 1); msk = msk(:,:,slice);
%             I_comb = maskon(cat(3, It1, It2));
%             It1 = I_comb(:,:,1);
%             It2 = I_comb(:,:,2);
%             It2 = imresize(It2, 0.5, "bilinear");
            imName = "slice" + string(slice);
    end
    % cropping the image
    [row, col] = ind2sub(size(It1), find(It1(:) > 0 & It2(:) > 0));
    
    crop = [min(row), max(row); min(col), max(col)];
%     % forcing image to have even dimensions
%     % this is ridiculous tbh
%     if any(mod([-1 1]*crop, 2))
%         b = mod([-1 1]*crop, 2)*[2;3];
%         switch b
%             case 2
%                 crop(1, 2) = crop(1, 2) + 1;
%             case 3
%                 crop(2, 2) = crop(2, 2) + 1;
%             case 5
%                 crop(:, 2) = crop(:, 2) + [1;1];
%         end
%     end
    It1 = It1(crop(1,1):crop(1, 2), crop(2,1):crop(2, 2));
    It2 = It2(crop(1,1):crop(1, 2), crop(2,1):crop(2, 2));
    msk = msk(crop(1,1):crop(1, 2), crop(2,1):crop(2, 2));
%     crop2 = [41, 41 + 89; 21, 21 + 83]; % crop 2: more crop
%     It1 = It1(crop2(1,1):crop2(1, 2), crop2(2,1):crop2(2, 2));
%     It2 = It2(crop2(1,1):crop2(1, 2), crop2(2,1):crop2(2, 2));
    images = struct();
    images.imName = imName;
    images.source = It2;
    hIt1 = It1;
    images.target = hIt1;
    
    [nx, ny] = size(It2);
    n = nx*ny;
    images.hX = reshape(hIt1, n, []);
    imDim = ceil([nx, ny] ./ patchSize);
    
    index = struct();
    images.crop = crop;
    images.slice = slice;
    if repType == "sample"
        [midx, grpLabel] = patchLabel([nx, ny], patchSize, repPoints);
        if numel(midx) == (imDim(1)-1)*imDim(2)
            imDim = imDim - [1 0];
        elseif numel(midx) == imDim(1)*(imDim(2)-1)
            imDim = imDim - [0 1];
        elseif numel(midx) == (imDim(1) - 1)*(imDim(2)-1)
            imDim = imDim - 1;
        end
        Y = double(reshape(It2, [], 1));
        
        It1 = reshape(hIt1, [], 1);
        It1 = It1(midx, :);
        It1 = reshape(It1, imDim(1), imDim(2), 1);
        images.rep = It1;
        hX = double(reshape(hIt1, [], 1));
        X = double(reshape(It1, [], 1));
        btidx = (hX ~= 0) & (Y ~= 0);
        if useFusedImage
            fY = imresize(It1, [nx ny], interType, "Antialiasing", false);
            fY = reshape(fY, size(Y, 1), []);
            Y = double([Y, fY]);
        end
        % binary and linear index for where at least one (BOTH is probably better) of the two images is
        % nonzero
        tidx = find(btidx);
        bigIndex = (1:n)'; % for any miscellaneous indexing needs
        
        
        bmidx = (1:n)' == grpLabel; % binary index of rep points
        gidx = find(~bmidx); % the (linear index) complement of midx
        mtidx = find(btidx & bmidx); % indices for truncated rep points
        sizevec = [reshape((repmat((1:ny)', [1, nx]))', [], 1), reshape(repmat((1:nx)', [1, ny]), [], 1)];
        sizevecColor = zeros(size(sizevec));
        sizevecColor(bmidx, :) = [reshape((repmat((1:imDim(2))', [1, imDim(1)]))', [], 1), reshape(repmat((1:imDim(1))', [1, imDim(2)]), [], 1)];
        
%         muT1 = 1/8; % normalization factors
%         muT2 = 1/8;
        muT1 = useSpatial*muX; muT2 = useSpatial*muY;
        % to turn off spatial info
        if muT1 == 0
            addX = zeros(nx*ny, 0);
        else
            addX = muT1.*(sizevecColor*diag(1./[imDim(2), imDim(1)]'));
        end
        
        if muT2 == 0
            addY = zeros(n, 0);
        else
            addY = muT2.*(sizevec*diag(1./[ny, nx]));
        end
        bkidx = logical(msk(:)); % size=[n, 1]
        btidxSmall = btidx(bmidx); % 1 for small points which are nz
        bkidxSmall = bkidx(bmidx); % 1 for small points which are mask
        bmskSmall = bkidxSmall(btidxSmall); % size=[size(X, 1), 1], 1 for rep pts which are mask
        target = [hX(bmidx & btidx, :), addX(bmidx & btidx, :)];
        rep = [Y(bmidx & btidx, :), addY(bmidx & btidx, :)];
        source = [Y(btidx, :), addY(btidx, :)];
        index.rep = bmidx;
        index.nz = btidx;
        index.mask = logical(msk(:));
        index.maskSmall = bmskSmall; % this is boolean
    elseif repType == "average"
        hX = double(reshape(hIt1, [], 1));
        Y = double(reshape(It2, [], 1));
        btidx = (hX ~= 0) & (Y ~= 0);
        It1 = imresize(hIt1, imDim, 'bilinear', 'Antialiasing', false);
        mskInter = imresize(double(msk), imDim, 'bilinear', 'Antialiasing', false);
        bkidxSmall = mskInter(:) >= 1;
        images.rep = It1;
        X = double(reshape(It1, [], 1));
        if useFusedImage
            fY = imresize(It1, [nx ny], interType, "Antialiasing", false);
            rfY = reshape(fY, size(Y, 1), []);
            Y = double([Y, rfY]);
        end
%         Z = imresize(It2, imDim, 'bilinear', 'Antialiasing', false);
        if useFusedImage
            Z = imresize(cat(3, It2, fY), imDim, 'bilinear', 'Antialiasing', false);
%             Z = imresize(cat(3, It2, fY), imDim, 'bicubic', 'Antialiasing', false);
        else
            Z = imresize(It2, imDim, 'bilinear', 'Antialiasing', false);
        end    
        Z = double(reshape(Z, prod(imDim), []));
        btidxSmall = (X ~= 0) & (sum(Z, 2) ~= 0);
        bmskSmall = bkidxSmall(btidxSmall);

        sizevec = [reshape((repmat((1:ny)', [1, nx]))', [], 1), reshape(repmat((1:nx)', [1, ny]), [], 1)];
        sizevecSmall = [reshape((repmat((1:imDim(2))', [1, imDim(1)]))', [], 1), reshape(repmat((1:imDim(1))', [1, imDim(2)]), [], 1)];
        
        muT1 = useSpatial*muX; muT2 = useSpatial*muY;
        % to turn off spatial info
        if muT1 == 0
            addX = zeros(prod(imDim), 0);
        else
            addX = reshape(sizevec, nx, ny, []);
            addX = imresize(addX, imDim, 'bilinear', 'Antialiasing', false);
%             addX = muT1.*reshape(addX, prod(imDim), [])*diag(1./[imDim(2), imDim(1)]');
            addX = muT1.*reshape(addX, prod(imDim), [])*diag(1./[ny nx]');

        end
        
        if muT2 == 0
            addY = zeros(n, 0);
        else
            addY = muT2.*(sizevec*diag(1./[ny, nx]));
        end
        target = [X(btidxSmall, :), addX(btidxSmall, :)];
        source = [Y(btidx, :), addY(btidx, :)];
        rep = [Z(btidxSmall, :), addX(btidxSmall, :)];
        index.rep = [];
        index.nz = btidx;
        index.mask = logical(msk(:));
        index.maskSmall = bmskSmall; % this is boolean
    end
end

end