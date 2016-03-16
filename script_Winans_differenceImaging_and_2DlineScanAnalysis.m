%%Script to determine mobile vesicle measurements for Winans et al. eLife
%%2016

%Written by A. Winans

%Attempted to different methods of background subtraction. First, used two 
%rounds of processing: imtophat filter (with diameter 1)
%followed by median filter. Second, used three rounds of processing: imtophat
%filter(with diameter 2), then median filter, then imtophat filter again
%(with diameter 1). This appears to greatly cut down on background.  


close all
clear all

readPath = '/Users/amy/Documents/Data/Nipkow/2014_02_21_Ftractin-mCherry&YFP-Syph/';

folders = {'cell1_18min_1', 'cell7_2', 'cell17_4min_1'};

window = 20;

load datMat


for j = 1
    
    files = dir(fullfile(readPath, folders{j}));
    names = arrayfun(@(x) x.name, files, 'UniformOutput', false);
    names = sort_im(names, 'img');
    
    namesRFP = sort_im(names, 'RFP');
    namesYFP = sort_im(names, 'YFP');
    
    imgStack = [];
    
    % read in all images
    for i = 1:length(namesRFP)
        
        img = imread(fullfile(readPath, folders{j}, namesRFP{i}));
        %     img = mat2gray(img, [0 2^12]);
        imgStack = cat(3, imgStack, img);
        
    end
    
    for i = 1:length(namesYFP)
        
        img = imread(fullfile(readPath, folders{j}, namesYFP{i}));
        %     img = mat2gray(img, [0 2^12]);
        imgStack = cat(3, imgStack, img);
        
    end
    
    
    % crop and mask
    imSum = sum(imgStack, 3);
    
    figure(9)
    imagesc(imSum)
    
 
        [crop, rect] = imcrop();

    
    
  crop2 = mat2gray(crop,[prctile(crop(:),5) prctile(crop(:),99.5)]);
    
    
    BW = minerrthreshold(crop2);
    
    BW2 = bwareaopen(BW, 1000);
    
    BW2 = bwmorph(BW2, 'close', 1);
    keyboard

        
        
        skel = bwmorph(BW2, 'skel', 'Inf');
        
        
        % Find single line path through neurite
        
        figure(10)
        imshow(skel)
        [x, y] = ginput(2);
        
        c1 = round(x(1));
        c2 = round(x(2));
        
        r1 = round(y(1));
        r2 = round(y(2));
               
        points{i} = [c1, c2; r1, r2];
        
        D1 = bwdistgeodesic(skel, c1, r1, 'quasi-euclidean');
        D2 = bwdistgeodesic(skel, c2, r2, 'quasi-euclidean');
        
        D = D1 + D2;
        D = round(D * 8) / 8;
        
        D(isnan(D)) = inf;
        skeleton_path = imregionalmin(D);
        
        [x, y] = find(skeleton_path); % single pixel path
                
        xcoord = x(1:window:end);
        ycoord = y(1:window:end);

    % Create windows in mask for each parameterized point
    
    cellMask = BW2;
    
    stats = regionprops(cellMask, 'Orientation');
    
    if abs(stats.Orientation) > 45
        [xcoord2, ind] = sort(xcoord);
        ycoord2 = ycoord(ind);
    else
        [ycoord2, ind] = sort(ycoord);
        xcoord2 = xcoord(ind);
    end
    
    windowCoord = [xcoord2, ycoord2];
    
    %apply filters to images
    
    imgStackF = [];
    int = [];
    for i = 1:10
        
        
        YFP = imgStack(:, :, i+1);
        
        YFP = imcrop(YFP, rect);
        
        YFPfilt = imtophat(YFP, strel('disk', 2));
        
        YFPfilt2 = medfilt2(YFPfilt);
        YFPfilt3 = imtophat(YFPfilt2, strel('disk', 1));
        
        axis image
        
        imgStackF = cat(3, imgStackF, YFPfilt3);
        
    end
    
    % take difference between images
    
    diffStack = [];
    close all
    for i = 1:9
        
        img1 = imgStackF(:, :, i);
        img2 = imgStackF(:, :, i+1);
        
        diff = img1 - img2;
        
        diff(diff < 0 ) = 0;
        
        diffStack = cat(3, diffStack, diff);
        
    end
    
    diffSum = sum(diffStack, 3);
    diffSum(diffSum  < 5 ) = 0;
    
    diffSum = uint16(diffSum);
    
        imwrite(diffSum, ['snapshot_syphdiff_t_' num2str(j) '.tif'], 'tif')
        
    % Create windows in mask for each parameterized point
    
    cellMask = BW2;
    
    stats = regionprops(cellMask, 'Orientation');
    
    if abs(stats.Orientation) > 45
        [xcoord2, ind] = sort(xcoord);
        ycoord2 = ycoord(ind);
    else
        [ycoord2, ind] = sort(ycoord);
        xcoord2 = xcoord(ind);
    end
    
    windowCoord = [xcoord2, ycoord2];
    [labelMask, boundaryOut] = getWindowLabelMap2(cellMask, windowCoord);
    
    RFP = imgStack(:, :, 1);
    
    RFP = imcrop(RFP, rect);
    
    imwrite(RFP, ['snapshot_act_' num2str(j) '.tif'], 'tif')
    
    RFP = mat2gray(RFP, [0 2^12]);
    RFP = background_subtract(RFP);
    actCrop = RFP;
    
            
    % Analysis per window
    
    win = max(labelMask(:));
    
    %actin window analysis
    
    actAll = [];
    
    for k = 1:win
        
        actWin = actCrop(labelMask == k);
        
        sorted = sort(actWin);
        highest = sorted(end-round(0.2*length(sorted)):end);
        medInt = median(double(highest));
        actAll = [actAll; medInt];
        
    end
    
    tujCrop = diffSum;
    
    tujAll = [];
    
    for k = 1:win
        
        tujWin = tujCrop(labelMask == k);
        
        tot = sum(tujWin);
        siz = length(tujWin);
        tujAll = [tujAll; [tot siz]];
      
    end
    
    figure(1)
    colormap('jet')
    
    subplot(3, 1, 1)
    imagesc(actCrop)
    hold on
    plot(y, x, 'r.', 'MarkerSize', 15);
    plot(c1, r1, 'g*', 'MarkerSize', 15)
    plot(c2, r2, 'g*', 'MarkerSize', 15)
    hold off
    axis image
    
    subplot(3, 1, 2)
    imagesc(diffSum)
    axis image
    set(gca, 'XTick', [])
    set(gca, 'YTick', [])
    
    subplot(3, 1, 3)
    imagesc(labelMask)
    set(gca, 'XTick', [])
    set(gca, 'YTick', [])
    axis image
    saveas(gcf, ['Snapshot_' num2str(j) '_Masking.jpeg'], 'jpeg')
    
    datMat(j).name = folders{j};
    datMat(j).rect = rect;
    datMat(j).BW = BW2;
    datMat(j).labelMask = labelMask;
    datMat(j).path = [x, y];
    datMat(j).dat = [actAll tujAll];
    
     save datMat datMat
    
end

% %%
%
% subplot(1, 3, 1)
% imagesc(YFP);
% axis image
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])
% subplot(1, 3, 2)
% imagesc(YFPfilt2)
% axis image
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])
% subplot(1, 3, 3)
% imagesc(sum(diffStack, 3));
% axis image
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])








