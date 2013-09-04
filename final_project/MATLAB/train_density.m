%% get the image names in the training folder
dir_name = 'dataset/train';
files = dir(fullfile([pwd '/' dir_name '/images'], '*.ppm'));
nimages = length(files);
fprintf('Total %d files.', nimages);

%% for each image, compute the superpixel features
nsp = 0;
features_all = [];
labels_all = [];
for i = 1 : nimages
    % load the image
    im = im2double(imread([dir_name '/images/' files(i).name]));
    grayim = rgb2gray(im);
    
    % load the image segmentation
    load([dir_name '/data/gt/' files(i).name '.mat']);
    
    % get the superpixel features
    spdata = get_spdata(im, grayim, imsegs);
    features = spdata.features;
    nsp = nsp + size(features, 1);
    
    features_all = [features_all; features];
    labels_all = [labels_all; imsegs.gvs];
    
    clc;
    fprintf('# superpixels: %d. Finished %s.\n', nsp, files(i).name);
end

%% then randomly select same-label and diff-label pairs
m = 1000;
% 1-1
ind = find(labels_all == 1);
ind = ind(randperm(length(ind)));
ind1 = ind(1 : m);
ind2 = ind(m + 1 : 2 * m);

features11 = abs(features_all(ind1, :) - features_all(ind2, :));

% 2-2
ind = find(labels_all == 2);
ind = ind(randperm(length(ind)));
ind1 = ind(1 : m);
ind2 = ind(m + 1 : 2 * m);

features22 = abs(features_all(ind1, :) - features_all(ind2, :));

% 3-3
ind = find(labels_all == 3);
ind = ind(randperm(length(ind)));
ind1 = ind(1 : m);
ind2 = ind(m + 1 : 2 * m);

features33 = abs(features_all(ind1, :) - features_all(ind2, :));

% 1-2
ind1 = find(labels_all == 1);
ind1 = ind1(randperm(length(ind1)));
ind1 = ind1(1 : m);

ind2 = find(labels_all == 2);
ind2 = ind2(randperm(length(ind2)));
ind2 = ind2(1 : m);

features12 = abs(features_all(ind1, :) - features_all(ind2, :));

% 1-3
ind1 = find(labels_all == 1);
ind1 = ind1(randperm(length(ind1)));
ind1 = ind1(1 : m);

ind2 = find(labels_all == 3);
ind2 = ind2(randperm(length(ind2)));
ind2 = ind2(1 : m);

features13 = abs(features_all(ind1, :) - features_all(ind2, :));

% 2-3
ind1 = find(labels_all == 2);
ind1 = ind1(randperm(length(ind1)));
ind1 = ind1(1 : m);

ind2 = find(labels_all == 3);
ind2 = ind2(randperm(length(ind2)));
ind2 = ind2(1 : m);

features23 = abs(features_all(ind1, :) - features_all(ind2, :));

% concatenate
features = [features11; features22; features33; features12; features13; features23];
labels = [ones(m, 1); ones(m, 1); ones(m, 1); -ones(m, 1); -ones(m, 1); -ones(m, 1)];

%% train the pairwise likelihood function
nfeatures = size(features, 2);

density = train_boosted_kde_2c(features, labels, [], 20);

for f = 1 : nfeatures
    ind = find(density(f).x > 0);
    density(f).x = density(f).x(ind)';
    density(f).log_ratio = density(f).log_ratio(ind)';
end

save('data/density.mat', 'density');