% runScripts
clear all;clc;close all;
numImages = 10;
matDir = 'Calib_Results_stereo.mat';
nr = 1440;
nc = 1920;
dirName = 'C:\Research\opencv\takePicture\images\depth\';
addpath(genpath('TOOLBOX_calib'));
codeMap = cell(1,15);


for depth = 1:75
    for lr = 0:1    
        getInvImage(dirName, numImages, lr, depth);
    end

    rectifyImage(matDir, dirName, numImages, depth);

    [codeL, minL, maxL, colorCodeL] = greyCode(nr, nc, numImages, 0, depth, dirName);

    fprintf('Code done for left...\n');
    [codeR, minR, maxR, colorCodeR] = greyCode(nr, nc, numImages, 1, depth, dirName);
    fprintf('Code done for right...\n');
    %idxL = (maxL-minL)/255 <= 1/255;
    %idxL = bwmorph(bwmorph(idxL, 'dilate', 3), 'erode', 3);
    idxL = (maxL-minL) > 4;
    bwmorph(idxL, 'clean');
    
    codeMap{depth}.codeL = colorCodeL;
    codeMap{depth}.idxL = idxL;
    
    idxR = (maxR-minR) > 4;
    bwmorph(idxR, 'clean');    
    
    codeMap{depth}.codeR = colorCodeR;
    codeMap{depth}.idxR = idxR;
    %idxR = (maxR-minR)/255 <= 1/255;
    %idxR = bwmorph(bwmorph(idxR, 'dilate', 3), 'erode', 3);
    
    codeL = codeL.*repmat(idxL, [1 1 size(codeL, 3)]);
    codeR = codeR.*repmat(idxR, [1 1 size(codeR, 3)]);
    
    
    for rr = 1 : nr
        if (sum(idxL(rr,:))  < 10)
            dispL(rr, :) = -(1:nc);
            dispR(rr, :) = -(1:nc);
        else
            scanL = codeL(rr,:, :); 
            scanL = reshape(scanL, size(codeL, 2), 1, numImages);
            scanR = codeR(rr, :,:); 
            

            [decL, decR]= convertGray(scanL, scanR);
            decL = repmat(decL, [1, size(scanR, 2)]);
            decR = repmat(decR, [size(scanL, 1),1]);
            dis = (decL - decR).^2;
%             scanR = reshape(scanR, 1, size(codeR, 2), numImages);
%             scanL =  repmat(scanL,[1 size(scanR, 2) 1]);
%             scanR = repmat(scanR, [size(scanL, 1) 1 1]);
%             
%             dis = sum(bitxor(scanL, scanR), 3);
            
            [~, tmpL] = min(dis, [], 2);
            [~, tmpR] = min(dis, [], 1);
            dispL(rr, :) = tmpL(:)' - (1:nc);
            dispR(rr, :) = tmpR(:)' - (1:nc);
            imagesc([dispL dispR]);colormap('jet');
            drawnow
        end
    end
%     dispL = medfilt2(dispL,[3 3]);
%     dispR = medfilt2(dispR,[3 3]);
    %dispL = dispL - repmat(1:nc, nr, 1);
    %dispR = dispR - repmat(1:nc, nr, 1);
    
%     imagesc([dispL dispR]); colormap('jet');
    save(sprintf('data/depth/disp_%02d.mat',depth), 'dispL', 'dispR');
    fprintf('Done....for depth %d', depth);
end
save('codeMapV1.mat', 'codeMap','-v7.3');
%% show the depthMap

% for i = 1%:10
%     file = sprintf('disp_%d.mat',i);
%     load(file);
%     %h = figure(i);
%     imagesc(abs([dispL - repmat(1:nc, nr, 1) dispR - repmat(1:nc, nr, 1)]), [125 325]);
%     %pause;
% %     saveas(h,sprintf('disp%d.png',i));
% end