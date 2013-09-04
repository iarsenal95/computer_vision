function maps = sp2regions(density, spdata, npart)

% cluster the superpixels into segments according to pairwise likelihoods
% given by density

features = spdata.features;

nsp = size(features, 1);
nfeatures = size(features, 2);
ndist = nsp * (nsp - 1) / 2;
nmaps = length(npart);

%% gets the matrix of pairwise sp likelihoods

% Y contains all pairwise likelihoods of being in same cluster
Y = zeros(ndist, 1);
for f = 1 : nfeatures
    tmpY = pdist(features(:, f), 'cityblock');
    % minus used because log ratio is of P(x|+)/P(x|-), but we want inverse    
    % NOW plus used for similarity
    Y = Y + get_kde_likelihood(density(f).log_ratio, density(f).x, tmpY);
end

Z = squareform(Y');
Z = Z .* (ones(nsp, nsp) - eye(nsp));

% convert Z to P(y1=y2 | |x1-x2|)
Z = log(1 ./ (1 + exp(-Z)));


%%

maps = zeros(nsp, nmaps);
for m = 1 : nmaps
    
    nclusters = npart(m);
    if nclusters > nsp
        nclusters = nsp;
    end
    
    % get random ordering of indices
    pind = randperm(nsp);
       
    % assign cluster ids to first nclusters clusters    
    for i = 1 : nclusters
        maps(pind(i), m) = i;
    end

    % assign remaining clusters
    % go forwards ...
    for i = nclusters + 1 : length(pind)
        % compute P(c(i) = k | x) for each k
        f = zeros(nclusters, 1);
        likelihoods = Z(pind(1 : i - 1), pind(i));            
        for k = 1 : nclusters
            sameinds = find(maps(pind(1 : i - 1), m) == k);         
            f(k) = mean(likelihoods(sameinds));
        end
       
        [~, k] = max(f);
        
        % assign cluster and add log(P(ci = k | x, c))
        maps(pind(i), m) = k;
    end
    
    %... and go backwards now using all cluster assignments
    for i = 1 : length(pind) % allow initial assignments to be reassigned
        % compute P(c(i) = k | x) for each k
        f = zeros(nclusters, 1);
        likelihoods = Z(:, pind(i));            
        for k = 1 : nclusters
            sameinds = find(maps(:, m) == k);          
            if ~isempty(sameinds)
                f(k) = mean(likelihoods(sameinds));
            end
        end 

        [~, k] = max(f);

        % assign cluster and add log(P(ci = k | x, c))
        maps(pind(i), m) = k;     
    end
    
    currmap = maps(:, m);
    nc = max(currmap);
    for k = 1 : nc
        ninds = length(find(currmap == k));
        while ninds == 0 && k < nc
            inds = find(currmap > k);
            currmap(inds) = currmap(inds) - 1;
            ninds = length(find(currmap == k));
            nc = max(currmap);
        end
    end
    maps(:, m) = currmap;

end