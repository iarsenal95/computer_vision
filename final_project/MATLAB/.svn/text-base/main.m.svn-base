clc;
clear all;
%% add path
addpath(genpath(pwd));

%% load the input image
imname = 'street07.ppm'; % TODO: modifiy the input image here
segname = 'tmpsp.ppm';
im = im2double(imread(imname));
% grayscale the image
grayim = rgb2gray(im);

%% segment the image
tic;

segcmd = 'segment/segment 0.8 100 100';
system([segcmd ' ' imname ' ' segname]);

%% pre-process the segmented image
imsegs = process_segments(segname);
imsegs.imname = imname;
% delete superpixel image
delete(segname);

%% then get superpixel data
disp('Computing superpixel features...');
spdata = get_spdata(im, grayim, imsegs);
disp('Done!');

%% find long straight edges
disp('Estimating horizon position...');
[lines spinfo] = find_long_edges(grayim, imsegs);

%% estimate the vanishing points
vpdata = estimate_vp(lines, size(grayim));
vpdata.lines = lines;
vpdata.spinfo = spinfo;

%% estimate the horizon position
vpdata.hpos = estimate_horizon(vpdata, size(grayim));
disp('Done!');

%% load the trained classifiers
load('data/classifier.mat');
load('data/density.mat');

%% group the superpixels into regions
disp('Grouping superpixels into regions...');
npart = [3 4 5 7 9 11 15 20 25];
nmaps = length(npart);
maps = sp2regions(density, spdata, npart);
disp('Done!');

%% now get the labels of superpixels(most time consuming)
disp('Labeling superpixels...');
labels = get_labels(im, imsegs, maps, spdata, vpdata, classifier);
disp('Done!');

%% write the VRLM file
disp('Writing VRLM model...');
write_vrlm_model(pwd, pwd, imsegs, labels, vpdata.hpos);
disp('Done!');

toc;