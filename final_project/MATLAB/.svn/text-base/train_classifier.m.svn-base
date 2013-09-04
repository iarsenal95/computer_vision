%% image names
dir_name = 'dataset/test';
files = dir(fullfile([pwd '/' dir_name '/images'], '*.ppm'));
nimages = length(files);
fprintf('Total %d files.\n', nimages);

load('data/density.mat');
npart = [3 4 5 7 9 11 15 20 25];
nreg = sum(npart);

%% for each image, extract region features
for i = 1 : nimages
    % load the image
    im = im2double(imread([dir_name '/images/' files(i).name]));
    grayim = rgb2gray(im);
    
    % load the image segmentation
    clear imsegs;
    load([dir_name '/data/gt/' files(i).name '.mat']);
    
    %% add the seginds field
    if ~isfield(imsegs, 'seginds')
        ns = imsegs.nseg;
        sinds = cell(ns, 1);
        for s = 1 : ns
            sinds{s} = find(imsegs.segimage == s);
        end
        imsegs.seginds = sinds;
    end
    
    %% get the superpixel features
    disp('Computing superpixel features...');
    spdata = get_spdata(im, grayim, imsegs);
    
    %% find long straight edges
    disp('Estimating vanishing points...');
    [lines spinfo] = find_long_edges(grayim, imsegs);

    %% estimate the vanishing points
    vpdata = estimate_vp(lines, size(grayim));
    vpdata.lines = lines;
    vpdata.spinfo = spinfo;

    %% estimate the horizon position
    vpdata.hpos = estimate_horizon(vpdata, size(grayim));
    
    %% group the superpixels to regions
    disp('Grouping superpixels...');
    maps = sp2regions(density, spdata, npart);
    nmaps = size(maps, 2);
    
    %% get region labels
    disp('Getting region labels...');
    rlabels = zeros(nreg, 1);
    k = 1;
    for m = 1 : nmaps
        curmap = maps(:, m);
        rlabels(k : k + npart(m) - 1) = get_region_labels(imsegs.gvs, curmap, npart(m));
        k = k + npart(m);
    end

    %% get region features
    disp('Computing region features...');
    rfeatures = zeros(nreg, 82);
    k = 1;
    for m = 1 : nmaps 
        curmap = maps(:, m);
        curnsegs = max(curmap(:));
        % get region features
        rfeatures(k : k + npart(m) - 1, :) = get_region_features(im, imsegs, curmap, (1 : curnsegs), spdata, vpdata);
        k = k + npart(m);
    end
    
    %% save to file
    disp('Saving...');
    save([dir_name '/data/rg/' files(i).name '_rg.mat'], 'rlabels', 'rfeatures', 'maps');
    
    %% print out
    clc;
    fprintf('%d / %d finished.\n', i, nimages);
end

%% now concatenate these data
nsamples = nreg * nimages;
data = zeros(nsamples, 82);
labels = zeros(nsamples, 1);
w = zeros(nsamples, 1);

cur = 0;
for i = 1 : nimages
    clear rfeatures;
    clear rlabels;
    load([dir_name '/data/rg/' files(i).name '_rg.mat']);
    data(cur + 1 : cur + nreg, :) = rfeatures;
    labels(cur + 1 : cur + nreg) = rlabels;
    w(cur + 1 : cur + nreg) = rfeatures(:, 41); % initial weights are proportional to regions' area
    cur = cur + nreg;
end

% normalize
w = w / sum(w);

%% train the classifier
disp('Training classifier...');
names = {'000', '090', 'sky', 'mix'};
ntrees = 20;
nnodes = 8;
classifier = train_boosted_dt_mc(data, [], names(labels)', ntrees, nnodes, 0, w, names);

save('data/classifier.mat', 'classifier');