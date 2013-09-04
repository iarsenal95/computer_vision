function labels = get_labels_test(imsegs, npart, afeatures, maps, vclassifier)

% we don't use horizontal classifiers

nsegs = imsegs.nseg;
nmaps = size(maps, 2);
nvclasses = length(vclassifier.names) - 1;
pV = zeros(nsegs, nvclasses);

% compute the probability of each block having each possible label
cur = 0;
for m = 1 : nmaps
    
    curmap = maps(:, m);
    curnsegs = npart(m);
    
    % get region features
    rfeatures = afeatures(cur + 1 : cur + curnsegs, :);
    nregions = size(rfeatures, 1);

    % get probability of vertical classes P(y|x) for each region
    tmpV = test_boosted_dt_mc(vclassifier, rfeatures);
    tmpV = 1 ./ (1 + exp(-tmpV));
    
    % normalize probabilities so that each is P(label|~mixed)P(~mixed)
    % normalize probabilities and sum over maps
    for r = 1 : nregions
        indices = find(curmap == r);  
        tmpV(r, 1 : end - 1) = tmpV(r, 1 : end - 1) / sum(tmpV(r, 1 : end - 1));
        tmpV(r, 1 : end - 1) = tmpV(r, 1 : end - 1) * (1 - tmpV(r, end));
        for c = 1 : nvclasses        
            pV(indices, c) = pV(indices, c) + tmpV(r, c);
        end                         
    end
    
    cur = cur + curnsegs;
end

% re-normalize weighted vote from classifiers
for s = 1 : size(pV, 1)
    pV(s, :) = pV(s, :) / sum(pV(s, :));
end

% get label for each block with confidence
labels = zeros(nsegs, 1);
for s = 1 : nsegs
    [~, c] = max(pV(s, :));
    labels(s) = c;  
end

end