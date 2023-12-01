function info =  paramHyperskinFun(startPath, arg)
    arguments
        startPath {mustBeFolder} 
        arg.faceNumber (1,1) {mustBePositive,mustBeInteger} = 3
        arg.poseNumber (1,1) {mustBePositive,mustBeInteger} = 3
        arg.sf (1,1) {mustBePositive} = 0.125
        arg.muX (1,1) {mustBeNonnegative} = 0
        arg.sigma (1,1) {mustBeNonnegative} = 0
        arg.numWorkers (1,1) {mustBePositive,mustBeInteger} = 1
        arg.keepVars {mustBeNumericOrLogical} = 1
        arg.nodims (1,1) {mustBePositive,mustBeInteger} = 5
        arg.K (1,1) {mustBePositive,mustBeInteger} = 10
        arg.modName (1,1) {mustBeValidModName(arg.modName)} = 3
    end

    if isstring(arg.modName) || ischar(arg.modName)
        [~, arg.modName] = ismember(upper(string(arg.modName)), ["MSI", "NIR", "RGB", "VIS", "FULL"]);
    end

    cd(char(startPath))
    addpath(genpath(char(startPath)))
    addpath(genpath('LE-fusion'))

    info = struct();
    
    [~, compName] = system('hostname'); % for identification purposes
    compName = string(compName(1:end-1));
    info.compName = compName;
    info.startPath = startPath;
    info.inputArgs = arg;
    
    % needs: faceNumber, poseNumber, modNum, sf, muX, numWorkers, keepVars
    
    
    %%
    dataSite = "https://www.comm.utoronto.ca/~pcng/data/ICASSP2024-SPGC/";
%     saveFolder = "/export/jemidih/"; % this is the directory that you want to save the files into, leave as empty string to save to current directory
    if isempty(startPath)
        startPath = "./";
    end
    saveFolder = startPath;
    modSite = ["Hyper-Skin(MSI,%20NIR)/", "Hyper-Skin(RGB,%20VIS)/"];
    msiFolder = ["MSI_CIE/", "NIR/"];
    rgbFolder = ["RGB_CIE/", "VIS/"];
    modFolder = [modSite(1)+ "train/" + msiFolder, modSite(2)+ "train/" + rgbFolder];
    
    
    
    % choose what face/pose/modality to care about
    % faceNumber = 3; % faces 3 through 51
    % poseNumber = 3; % 6 poses, 3 neutrals and 3 smiles, front left right, in abc order
    % modNum = 4; % four modality types in abc order in terms of (information content), multispectral (high), near IR (low), rgb (low), vis (high)
    
    % sprintf("p%03d", n)
    recycle('off')
    poseList = reshape(["neutral_", "smile_"] + ["front", "left", "right"]', 1, []);
    fileName = sprintf("p%03d_%s.", arg.faceNumber, poseList(arg.poseNumber));
    if arg.modName ~=5
        dataDir = modFolder(arg.modName);
        dataUrl = dataSite + dataDir;
    end
    options = weboptions('RequestMethod','get', 'ContentType', 'image');
    
    % A = [];
    
    % A = load(websave(dataUrl(2), dataUrl(1)+dataUrl(2), options));
    % make something to hold the data that has the same file structure, then do
    % like if ~exist, load(websave...)
    
    % set up file structure
    for i = 1:numel(modFolder)
        if ~exist(saveFolder + modFolder(i), 'dir')
            mkdir(saveFolder + modFolder(i))
        end
    end
    
    % % the rgb files are just jpegs
    % if modNum == 3
    %     fileName = replace(fileName, ".mat", ".jpg");
    % end
    
    % check if file exists
    if arg.modName == 5
        if exist(saveFolder + modFolder(4) + fileName + "mat", 'file')
            A = getfield(load(saveFolder + modFolder(4) + fileName + "mat"), 'cube');
        else
            A = getfield(load(websave(saveFolder + modFolder(4) + fileName + "mat", dataSite + modFolder(4) + fileName + "mat", options)), 'cube');
        end

        if exist(saveFolder + modFolder(2) + fileName + "mat", 'file')
            Aplus = getfield(load(saveFolder + modFolder(2) + fileName + "mat"), 'cube');
            A = cat(3, A, Aplus(:,:,2:end));
        else
            Aplus = getfield(load(websave(saveFolder + modFolder(2) + fileName + "mat", dataSite + modFolder(2) + fileName + "mat", options)), 'cube');
            A = cat(3, A, Aplus(:,:,2:end));
        end
    elseif arg.modName ~= 3
        if exist(saveFolder + dataDir + fileName + "mat", 'file')
            A = getfield(load(saveFolder + dataDir + fileName + "mat"), 'cube');
        else
            A = getfield(load(websave(saveFolder + dataDir + fileName + "mat", dataUrl + fileName + "mat", options)), 'cube');
        end
    else
        if exist(saveFolder + dataDir + fileName + "jpg", 'file')
            A = imread(saveFolder + dataDir + fileName + "jpg");
        else
            A = imread(websave(saveFolder + dataDir + fileName + "jpg", dataUrl + fileName + "jpg", options));
        end
    end
    
    % A is 1024x1024x31 (or 3 for RGB)
    if isinteger(A)
    A = double(A)./(pow2(8*(whos('A').bytes)/numel(A)) - 1);
    end
    % overengineered way to convert to a double from an int
    
    % if we want to resize by 1/2 (or some other 0 < scale < 1) in both length and width:
    % sf = 0.125;
    A = imresize(A, arg.sf);
    
    %% make and apply a mask to keep only the skin/brighter parts of the image, based on RGB data
    
    if exist(saveFolder + modFolder(3) + fileName + "jpg", 'file')
        msk = imread(saveFolder + modFolder(3) + fileName + "jpg");
    else
        msk = imread(websave(saveFolder + modFolder(3) + fileName + "jpg", dataSite + modFolder(3) + fileName + "jpg", options));
    end
    
    
    if isinteger(msk)
    msk = double(msk)./(pow2(8*(whos('msk').bytes)/numel(msk)) - 1);
    end
    
    msk = imresize(msk, size(A, 1:2));
    % msk = sum(msk, 3) > median(sum(msk, 3), 'all'); % hoping that most of the image is dark
    msk = filtering_b_channel(msk);
    
    %%
    [nx, ny] = size(A, [1 2]);
    sizevec = [reshape((repmat((1:ny)', [1, nx]))', [], 1), reshape(repmat((1:nx)', [1, ny]), [], 1)];
    % muX = 1/2;
    
    if arg.muX == 0
        addX = zeros(prod([nx, ny]), 0);
    else
        addX = arg.muX.*sizevec*diag(1./[ny, nx]');
    end
    
    
    %%
    X = reshape(A, [], size(A, 3)); 
    if 2*nnz(sum(X, 2)) <= prod(size(A, 1:2))
        X = sparse(X);
    end
    mskX = X(msk(:), :);
    mapX.X = [mskX, addX(msk(:), :)];
    % nodims = 40;
    % K = 100;
    % sigma = 0;
    [mapX.vec, mapX.val, ~, mapX.sigma, ~, mapX.aff] = lapEig(mskX, arg.nodims, arg.K, 'sigma', arg.sigma); % embed
    
    % preX = zeros(size([X, addX])); % preallocate
    preX = spalloc(size([X, addX], 1), size([X, addX], 2), nnz(msk)*size([X, addX], 3));
    t_start = tic;
    preX(msk(:), :) = LEPre_init(mapX.vec, mapX, 'numWorkers', arg.numWorkers); % then preimage, 
    t_end = toc(t_start);
    fprintf('Completed in %.2f seconds. \n', t_end)
    info.preimageTime = t_end;
    preX = preX(:, 1:size(A, 3)); % get rid of recovered spatial coords, if any
    
    preA = reshape(full(preX), size(A));
    
    ssimX = ssim(preA, A.*msk);
    psnrX = psnr(preA, A.*msk);
    info.ssim = ssimX;
    info.psnr = psnrX;
    
    
    badMask = any(preA > 1, 3) | any(preA < 0, 3) | any(isnan(preA), 3);
    badMask = ~badMask(:);
    
    fprintf("Structural similarity index: %0.4f \n", ssimX)
    fprintf("Peak signal-to-noise ratio: %2.3f \n", psnrX)
    fprintf("%d out of %d pixels in range. \n", nnz(badMask & msk(:)), nnz(msk(:)))
    info.pixelsInRange = sprintf("%d out of %d pixels in range", nnz(badMask & msk(:)), nnz(msk(:)));
    
    
    saveFile = sprintf("%dTest_Face%d_Pose%d_mod%d_"+string(datetime('now', 'Format', 'yy-MM-dd-HHmm')),...
        nx, arg.faceNumber, arg.poseNumber, arg.modName);
    
    
    if arg.keepVars
    info.saveFile = saveFile;
    save(saveFile, '-regexp', '^((?!map).)*$'),
    save("map_"+saveFile, '-regexp', '(map)\w+')
    else
        info.saveFile = "";
    end
    fprintf("Done! \n")
    


end

function mustBeValidModName(a)
    if numel(a) ~= 1
        eidType = 'mustBeValidModName:need1Args';
        msgType = 'Input must have 1 element. ';
        error(eidType,msgType)
    end
    if isstring(a) && ischar(a)
        if ~all(ismember(upper(string(a)), ["MSI", "NIR", "RGB", "VIS", "FULL"]))
            eidType = 'mustBeValidModName:badName';
            msgType = 'Input must be one of the 3-letter codes. ';
            error(eidType,msgType)
        end
    elseif isnumeric(a) && (~isequal(a, floor(a)) || any(a > 5) || any(a < 1))
            eidType = 'mustBeValidModName:badNum';
            msgType = 'Input must be either 1, 2, 3, 4, or 5. ';
            error(eidType,msgType)
    end
end