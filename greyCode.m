function [codeMap, minI, maxI, codeColor] = greyCode(nr,nc, numImages,lr, depth, dirName)

%dirName = 'C:\Research\opencv\takePicture\images\depth\';
dirFile = [dirName sprintf('%d\\',depth)];
if (lr)
    baseName = sprintf('R_rect_');
    revName = sprintf('R_rect_inv_');
else
    baseName = sprintf('L_rect_');
    revName = sprintf('L_rect_inv_');
end

codeColor = zeros(nr,nc);
codeMap = zeros(nr,nc, numImages);

minI = Inf;
maxI = 0;

for i = 1:numImages
    I1 = double(imread([dirFile baseName sprintf('%02d.png',i)]));
    I2 =  double(imread([dirFile revName sprintf('%02d.png',i)]));
    
    minI = min(minI, I1);
    minI = min(minI, I2);
    maxI = max(maxI, I1);
    maxI = max(maxI, I2);
    
    
    Itmp = I1 > I2;
    
     codeColor    = codeColor + Itmp * 2^(i); %numImages-
%     imagesc(Itmp);
%     pause;
    codeMap(:,:,i) = Itmp;
end