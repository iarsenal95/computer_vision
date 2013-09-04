function spdata = get_spdata(im, grayim, imsegs)

% load the pre-computed doog filters
load('data/doog_filters.mat');
nang = size(doog_filters, 3);

% convolve the image with each orientation
orientim = zeros(size(im, 1), size(im, 2), nang);
for i = 1 : nang
    orientim(:, :, i) = abs(conv2(grayim, doog_filters(:, :, i), 'same') - grayim);
end

% compute the features in each superpixel
nseg = imsegs.nseg;
nfeature = nang + 3 + 6 + 2;
features = zeros(nseg, nfeature);

% 1 - 12: mean doog response
maar = zeros(nseg, nang);
for i = 1 : nang
    maar(:, i) = get_spmeans(imsegs, orientim(:, :, i));
end

features(:, 1 : nang) = maar;

% 13: mean of mean doog response in all orientations
features(:, nang + 1) = mean(maar, 2);

% 14: id of max response of all orientations
[maxval, features(:, nang + 2)] = max(maar, [], 2);

% 15: (max - median) of maar
features(:, nang + 3) = maxval - median(maar, 2);

spdata.orientation = maar;
spdata.edginess = mean(orientim, 3);

% 16 - 18: means of r, g, b
rgb = zeros(nseg, 3);
% r, g, b channels
for i = 1 : 3
    rgb(:, i) = get_spmeans(imsegs, im(:, :, i));
end
features(:, nang + (4 : 6)) = rgb;

% 19 - 21: means of h, s, v
features(:, nang + (7 : 9)) = rgb2hsv(rgb);

% 23: normalized x and y
[h, w, ~] = size(im);
yim = 1 - repmat(((0 : h - 1) / (h - 1))', 1, w);
xim = repmat(((0 : w - 1)/(w - 1)), h, 1);
features(:, nang + 10) = get_spmeans(imsegs, yim);
features(:, nang + 11) = get_spmeans(imsegs, xim);

spdata.features = features;

% hsv
[spdata.hue, spdata.sat, ~] = rgb2hsv(im);

end