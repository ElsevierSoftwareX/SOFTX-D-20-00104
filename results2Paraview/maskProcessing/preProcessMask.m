function mask_out = preProcessMask(mask_in)


%Dilate
SE = strel('sphere',2);
mask_in = imdilate(mask_in, SE);
%Fill holes
mask_in = imfill(mask_in,'holes');
%Erode
mask_in = imerode(mask_in, SE); clear SE;
%Keep the largest mask objects
conn = 3^ndims(mask_in)-1;
[mask_in, ~] = getLargestCc(mask_in,conn,4);
%Add padding
padSize = 3;
mask_in=padarray(mask_in,[padSize padSize padSize]);
%Smooth out the mask surface
[mask_in,~] = moisturize (mask_in);
%Remove the previous padding
mask_out = mask_in(padSize+1:end-padSize,padSize+1:end-padSize,padSize+1:end-padSize);