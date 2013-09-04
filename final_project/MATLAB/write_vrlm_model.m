function write_vrlm_model(imdir, vrmldir, imseg, labels, hpos)

fn = imseg.imname;
image = im2double(imread([imdir '/' fn]));    

% horizon position
labels.hy = 1 - hpos;

% labels
vlabels = labels.vert_labels;
hy = labels.hy;

% eliminate the borders
bn = strtok(fn, '.');            
imseg.segimage = imseg.segimage(26:end-25, 26:end-25);
image = image(26:end-25, 26:end-25, :);

[gplanes, vplanes, gmap,vmap] = APPlabels2planes(vlabels, [], hy, imseg.segimage, image);

use_fancy_transparency = 1;
if use_fancy_transparency
    vmap = conv2(double(vmap), fspecial('gaussian', 7, 2), 'same');
    vim = image;
    vim(:, :, 4) = vmap; % add alpha channel
else
    vim = image;
    for b = 1:size(image, 3)
        vim(:, :, b) = image(:, :, b).*vmap;
    end
end

gim = image;    
for b = 1:size(image, 3)
    gim(:, :, b) = gim(:, :, b) .* gmap;
end    

[points3D, tpoints2D, faces] = APPplanes2faces(gplanes, vplanes, [size(image, 1) size(image, 2)], hy);   
faces2vrml(vrmldir, bn, points3D, tpoints2D, faces, gim, vim);
    
end