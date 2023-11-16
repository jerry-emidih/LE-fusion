% start from the folder that has the code in it
addpath(genpath('LE-fusion'))


%%
dataSite = "https://www.comm.utoronto.ca/~pcng/data/ICASSP2024-SPGC/";
saveFolder = ""; % this is the directory that you want to save the files into, leave as empty string to save to current directory
modSite = ["Hyper-Skin(MSI,%20NIR)/", "Hyper-Skin(RGB,%20VIS)/"];
msiFolder = ["MSI_CIE/", "NIR/"];
rgbFolder = ["RGB_CIE/", "VIS/"];
modFolder = [modSite(1)+ "train/" + msiFolder, modSite(2)+ "train/" + rgbFolder];



% choose what face/pose/modality to care about
faceNumber = 13; % faces 3 through 51
poseNumber = 3; % 6 poses, 3 neutrals and 3 smiles, front left right, in abc order
modNum = 3; % four modality types in abc order in terms of (information content), multispectral (high), near IR (low), rgb (low), vis (high)

% sprintf("p%03d", n)
poseList = reshape(["neutral_", "smile_"] + ["front", "left", "right"]', 1, []);
dataDir = modFolder(modNum);
fileName = sprintf("p%03d_%s.mat", faceNumber, poseList(poseNumber));
dataUrl = dataSite + dataDir;
options = weboptions('RequestMethod','get', 'ContentType', 'image');

% A = [];

% A = load(websave(dataUrl(2), dataUrl(1)+dataUrl(2), options));
% make something to hold the data that has the same file structure, then do
% like if ~exist, load(websave...)

% set up file structure
for i = 1:numel(modFolder)
    if ~exist(modFolder(i), 'dir')
        mkdir(modFolder(i))
    end
end

% the rgb files are just jpegs
if modNum == 3
    fileName = replace(fileName, ".mat", ".jpg");
end

% check if file exists
if modNum ~= 3
    if exist(dataDir + fileName, 'file')
        A = load(dataDir + fileName);
    else
        A = load(websave(dataDir + fileName, dataUrl + fileName, options));
    end
    A = A.cube;
else
    if exist(dataDir + fileName, 'file')
        A = imread(dataDir + fileName);
    else
        A = imread(websave(dataDir + fileName, dataUrl + fileName, options));
    end
end

% A is 1024x1024x31 (or 3 for RGB)

if isinteger(A)
    A = double(A)./(pow2(8*(whos('A').bytes)/numel(A)) - 1);
end
% overengineered way to convert to a double from an int


% if we want to resize by 1/2 (or some other 0 < scale < 1) in both length and width:
sf = 1/16; % scaling factor
A = imresize(A, sf);

X = reshape(A, [], size(A, 3)); 
% each pixel of A is a row of X, in column-major order

[mapX.vec, mapX.val, ~, mapX.sigma, ~, mapX.aff] = lapEig(X, 15, 50, 'sigma', 0);

%% Notes on mapX:
% .vec is the matrix of embedded coordinates
% .val is a vector of eigenvalues
%
% Here, we want the first 15 nonzero eigenvalues, so if 0 has a high
% multiplicity, the code will do multiple eigensolves until we get what we
% want (or it breaks)
%
% mapX.sigma is the bandwidth parameter sigma used in the weight function.
% in normLap, I have it so that with input 0, it sets sigma to the 80th
% percentile of the measured neighbor distances




