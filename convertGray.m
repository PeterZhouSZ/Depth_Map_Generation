function [decL, decR]= convertGray(grayL, grayR)
    nbits = 10;
    grayTemplate = generateGray(nbits);
    for i = 1:nbits
        grayTemplate(:,i) = grayTemplate(:,i)*2^(nbits - i);
        tempL(1,:,i) = grayL(:,:,i)*2^(nbits - i);
        tempR(1,:,i) = grayR(:,:,i)*2^(nbits - i);
    end
    grayTemplate = sum(grayTemplate, 2);
    tempL = sum(tempL,3);
    tempR = sum(tempR,3);
    for j = 1:length(tempL)
        decL(j,1) = find(grayTemplate == tempL(j));
        decR(1,j) = find(grayTemplate == tempR(j));
    end
end

function Gnew = generateGray(nbits)
% nbits = 10;

G = [ 0 0; 0 1; 1 1; 1 0];
for kk=3:nbits
    G = [ zeros(size(G, 1), 1) G; ones(size(G, 1), 1) G(end:-1:1, :)];
end

siz  = [ 1920 1080];
Gnew = G;%G(1:siz(1), :);
end
