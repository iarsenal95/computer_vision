%% load base names
dir_name = 'dataset/test';
files = dir(fullfile([pwd '/' dir_name '/images'], '*.ppm'));
nimages = length(files);
fprintf('Total %d files.\n', nimages);

load('data/classifier.mat');
npart = [3 4 5 7 9 11 15 20 25];

%% load ground truth and region features, and get the superpixel labels
nlabels = zeros(3, 1);
nerrors = zeros(3, 1);

for i = 1 : nimages
    clear imsegs;
    load([dir_name '/data/gt/' files(i).name '.mat']);
    clear rfeatures;
    load([dir_name '/data/rg/' files(i).name '_rg.mat']);
    
    labels = get_labels_test(imsegs, npart, rfeatures, maps, classifier);
    
    error = labels - imsegs.gvs;
    
    for j = 1 : 3
        ind = find(imsegs.gvs == j);
    	nlabels(j) = nlabels(j) + length(ind);
    	nerrors(j) = nerrors(j) + nnz(error(ind));
    end
    disp(i);
end

fprintf('Total accuracy: %f.\n', 1 - sum(nerrors) / sum(nlabels));
fprintf('Ground accuracy: %f.\n', 1 - nerrors(1) / nlabels(1));
fprintf('Vertical accuracy: %f.\n', 1 - nerrors(2) / nlabels(2));
fprintf('Sky accuracy: %f.\n', 1 - nerrors(3) / nlabels(3));