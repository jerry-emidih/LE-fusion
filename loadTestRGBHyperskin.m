

dataSite = "https://www.comm.utoronto.ca/~pcng/data/ICASSP2024-SPGC/";
startPath = [];
if isempty(startPath)
    startPath = "./";
end
saveFolder = startPath;
modSite = ["Hyper-Skin(MSI,%20NIR)/", "Hyper-Skin(RGB,%20VIS)/"];
msiFolder = ["MSI_CIE/", "NIR/"];
rgbFolder = ["RGB_CIE/", "VIS/"];
modFolder = [modSite(1)+ "test/" + msiFolder, modSite(2)+ "test/" + rgbFolder];


% n = 50; % number of training points
% faceIdx = [3:15, 17:29, 31:34, 37:51];
faceIdx = [1, 2, 3, 12, 16, 30, 35, 36, 39];
% faceIdx = [31:34, 37:51]';
[rows, cols] = ndgrid(1:6, faceIdx);
runID = [rows(:), cols(:)];

poseList = reshape(["neutral_", "smile_"] + ["front", "left", "right"]', 1, []);
options = weboptions('RequestMethod','get', 'ContentType', 'image');




modNum = 3; % rgb images
for ii = 1:size(runID, 1)
           
        recycle('off')
        poseList = reshape(["neutral_", "smile_"] + ["front", "left", "right"]', 1, []);
        dataDir = modFolder(modNum);
        fileName = sprintf("p%03d_%s.", runID(ii, 2), poseList(runID(ii, 1)));
        dataUrl = dataSite + dataDir;
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
        if modNum ~= 3
            if exist(saveFolder + dataDir + fileName + "mat", 'file')
                A = load(saveFolder + dataDir + fileName + "mat");
            else
                A = load(websave(saveFolder + dataDir + fileName + "mat", dataUrl + fileName + "mat", options));
            end
            A = A.cube;
        else
            if exist(saveFolder + dataDir + fileName + "jpg", 'file')
                A = imread(saveFolder + dataDir + fileName + "jpg");
            else
                A = imread(websave(saveFolder + dataDir + fileName + "jpg", dataUrl + fileName + "jpg", options));
            end
        end

end