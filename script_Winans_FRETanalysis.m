%% Sample script to describe image processing of FRET images. For Winans et
%% al., eLife 2016

% Written by A. Winans

%% Briefly,
%% images are pre-cropped, read into program, subject to a flat, median
%% background subtraction, and then filtered through a gaussian filter
%% before taking the ratio of FRET/CFP images.


close all
clear all

window = 5;

load datMatCdc
load datMatCdc2

readPath = '/Users/amy/Documents/MATLAB/scripts/actinwave/Rac1 sensor/dat082814_ixmicro/ss_Cdc42/';

folders = dir(readPath);

foldersName = {folders.name};

foldersName = sort_im(foldersName, 'ss');

for i = 15:length(foldersName)
    
    files = dir(fullfile(readPath, foldersName{i}));
    
    names = {files.name};
    
    FRETnames = sort_im(names, 'FRET');
    CFPnames = sort_im(names, 'CFP');
    RFPnames = sort_im(names, 'Texas');
     
    imFRET = double(imread(fullfile(readPath,foldersName{i}, FRETnames{1})));
    imCFP = double(imread(fullfile(readPath, foldersName{i},CFPnames{1})));
    imRFP = double(imread(fullfile(readPath, foldersName{i},RFPnames{1})));
    


bgMask0 = datMatCdc(i).BW;

bgMask0 = ~bgMask0;
bgMask = imerode(bgMask0, strel('disk', 3));


bgMaskedFRET = bgMask.*imFRET;
bgMaskedCFP = bgMask.*imCFP;
bgMaskedRFP = bgMask.*imRFP;


bgMaskedFRET(bgMaskedFRET == 0) = [];
bgMaskedCFP(bgMaskedCFP == 0) = [];
bgMaskedRFP(bgMaskedRFP == 0) = [];

bgFRET = median(bgMaskedFRET(:));
bgCFP = median(bgMaskedCFP(:));
bgRFP = median(bgMaskedRFP(:));

imFRET = imFRET - bgFRET;
imCFP = imCFP - bgCFP;
imRFP = imRFP - bgRFP;


   
BW = datMatCdc(i).BW;
    BW2 = BW;

    imFRETF = imFRET;
    imCFPF = imCFP;
    
    imFRETF(~BW)=nan;
    imFRETF=ndnanfilter(imFRETF,fspecial('gaussian',11,3),11);
    imFRETF(~BW)=nan;
    imCFPF(~BW)=nan;
    imCFPF=ndnanfilter(imCFPF,fspecial('gaussian',11,3),11);
    imCFPF(~BW)=nan;
    
    imRatio=imFRET./imCFP;
    imRatioF = imFRETF./imCFPF;
    
    figure(1)
    subplot(1, 2, 1)
    imagesc(imRFP, [0 3000])
    colormap('gray')
    axis image
    axis off
    
    subplot(1, 2, 2)
    imagesc(BW2.*imRatioF, [0 2])
    colormap('jet')
    axis image
    axis off
    
    saveas(gcf, ['Im5_' foldersName{i} '_Sensor.jpeg'], 'jpeg')
