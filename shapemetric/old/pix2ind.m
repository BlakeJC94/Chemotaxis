function ind = pix2ind(pix, im)
% converts pix (1x2 row,column) into index
ind = (pix(2)-1)*size(im,2) + pix(1);
end