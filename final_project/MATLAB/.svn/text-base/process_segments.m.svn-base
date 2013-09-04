function imsegs = process_segments(segname)

% the sturcture of segments
imsegs = struct('imsize', [0 0], 'segimage', [], 'nseg', 0, 'adjmat', []);

% load the segmented image
im = double(imread(segname));
imsegs.imsize = size(im);
imsegs.imsize = imsegs.imsize(1 : 2);

% 3 uint8 -> uint32 -> double
im = im(:, :, 1) + im(:, :, 2) * 256 + im(:, :, 3) * 256 * 256;

% label the segments
[gid, gn] = grp2idx(im(:));
imsegs.segimage = uint16(reshape(gid, imsegs.imsize));
imsegs.nseg = length(gn);

% compute the adjacency matrix according to the labels
nseg = imsegs.nseg;
segimage = imsegs.segimage;
h = imsegs.imsize(1);
adjmat = eye([nseg nseg], 'uint8');

% adj-diff along x and y axis
dx = uint8(segimage ~= segimage(:, [2 : end end]));
dy = uint8(segimage ~= segimage([2 : end end], :));

% vertical adj
ind1 = find(dy);
ind2 = ind1 + 1;
s1 = segimage(ind1);
s2 = segimage(ind2);
adjmat(s1 + nseg * (s2 - 1)) = 1;
adjmat(s2 + nseg * (s1 - 1)) = 1;

% horizontal adj
ind3 = find(dx);
ind4 = ind3 + h;
s3 = segimage(ind3);
s4 = segimage(ind4);
adjmat(s3 + nseg * (s4 - 1)) = 1;
adjmat(s4 + nseg * (s3 - 1)) = 1;

stats = regionprops(segimage, 'Area');
imsegs.npixels = vertcat(stats(:).Area);
imsegs.adjmat = logical(adjmat);

ns = imsegs.nseg;
sinds = cell(ns, 1);
for s = 1 : ns
    sinds{s} = find(imsegs.segimage == s);
end
imsegs.seginds = sinds;

end