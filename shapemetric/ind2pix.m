function pix = ind2pix(ind, im)
% converts pix (1x2 row,column) into index
pix = [mod(ind, size(im,2)), ...
    (ind - mod(ind, size(im,2)))/size(im,2) + 1];
end