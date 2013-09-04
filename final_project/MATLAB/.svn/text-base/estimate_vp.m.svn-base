function vpdata = estimate_vp(lines, imsize)

nlines = size(lines, 1);

x1 = [lines(:, [1 3]) ones(size(lines, 1), 1)];
x2 = [lines(:, [2 4]) ones(size(lines, 1), 1)];

% get plane normals for line segments
l = cross(x1, x2);
l = l ./ repmat(sqrt(sum(l .^ 2, 2)), 1, 3);

% make theta range from 0 to pi (instead of -pi to pi)
nbins = 60;
theta = mod(lines(:, 5), pi);

% get histogram of thetas
binwidth = pi / nbins;
bincenters = (binwidth / 2) : binwidth : (pi - binwidth / 2);
hist_theta = hist(theta, bincenters);

% smooth histogram
for b = 1 : nbins
    hist_theta(b) = sum(hist_theta(mod(((b - 1) : (b + 1)) - 1, nbins) + 1) .* [0.25 0.5 0.25]);
end

% compute curvature of histogram
C = zeros(1, nbins);
for b = 1:nbins
    C(b) = hist_theta(b) - mean(hist_theta(mod(((b - 4) : (b + 4)) - 1, nbins) + 1));
end


% find dominant peaks surrounded by zero crossings
zc_pos = find((C > 0) & ([C(end) C(1 : end - 1)] < 0));
zc_neg = find((C < 0) & ([C(end) C(1 : end - 1)] > 0));

ngroups = length(zc_pos);
groups = cell(1, ngroups);
bc = bincenters + pi / nbins / 2;
if zc_neg(1) < zc_pos(1)
    i1 = round((zc_pos(end) + zc_neg(end)) / 2);
    i2 = round((zc_neg(1) + zc_pos(1)) / 2);
    groups{1} = find((theta > bc(i1)) | (theta < bc(i2)));
    for i = 2 : ngroups
        i1 = round((zc_pos(i - 1) + zc_neg(i - 1)) / 2);
        i2 = round((zc_neg(i) + zc_pos(i)) / 2);     
        groups{i} = find((theta > bc(i1)) & (theta < bc(i2)));
    end

else
    for i = 1 : ngroups
        i1 = ceil((zc_pos(i) + zc_neg(max(i - 1, 1))) / 2);
        i2 = ceil((zc_neg(i) + zc_pos(min(i + 1, ngroups))) / 2);
        groups{i} = find((theta > bc(i1)) & (theta < bc(i2)));       
    end
end

% remove groups with few segments 
remove = [];
thresh = max(0.05 * nlines, 5);
for i = 1 : ngroups
    if length(groups{i}) < thresh
        remove(end + 1) = i;
    end
end
groups(remove) = [];
ngroups = length(groups);


% initialize EM
p = zeros(nlines, ngroups); % p = P(v | l) = P(l | v)P(v) / P(l)
for i = 1 : ngroups
    p(groups{i}, i) = 1;
end

A = l;
v = zeros(3, ngroups);
sigma = ones(1, ngroups); % variance
for i = 1 : ngroups
    normp = p(:, i) / sum(p(:, i));
    W = diag(normp);
    [eigV, lambda] = eig(A' * (W') * W * A);
    [~, smallest] = min(diag(lambda));
    v(:, i) = eigV(:, smallest);
    sp = sort(normp, 'descend');
    sp = sum(sp(1:min(length(sp), 2)));
    sigma(i) = normp' * (l * v(:, i)) .^ 2 / (1 - sum(sp));
end

% create outlier groups
nadd = 3;
ngroups = ngroups + nadd;
sigma(end + (1 : nadd)) = 0.2;
tmpv = [0 0 1]';
v(:, end + 1) = tmpv / sqrt(sum(tmpv .^ 2));
tmpv = [1 0 1]';
v(:, end + 1) = tmpv / sqrt(sum(tmpv .^ 2));
tmpv = [-1 0 1]';
v(:, end + 1) = tmpv / sqrt(sum(tmpv .^ 2));
p(:, end + (1 : nadd)) = 0;
pv = ones(1, ngroups);

oldv = v;

% do EM
for iter = 1:15
    pv = pv / sum(pv);
    S = repmat(sigma, nlines, 1) + 1E-10;
    plv = exp(-(l * v) .^ 2 ./ S / 2) ./ sqrt(S) + 1E-10;        
    p = plv .* repmat(pv, nlines, 1);    
    p = p ./ repmat(sum(p, 2), 1, ngroups);  
    pv = sum(p, 1);
    for i = 1 : ngroups
        normp = p(:, i) / sum(p(:, i));
        W = diag(normp);
        [eigV, lambda] = eig(A' * (W') * W * A);
        [~, smallest] = min(diag(lambda));
        v(:,i) = eigV(:, smallest);            
        sp = sort(normp, 'descend');
        sp = sum(sp(1 : min(length(sp), 2)));
        if (1 - sum(sp)) > 0
            sigma(i) = normp' * (l * v(:, i)) .^ 2 / (1 - sum(sp)); 
        else
            sigma(i) = Inf;                       
        end
    end        

    % remove duplicate groups
    remove =[];
    for i = 1 : ngroups        
        for j = i + 1 : ngroups
            if (v(:, i)' * v(:, j) > 0.995) || (pv(i) * nlines <= 3) || (sigma(i) > 10) 
                remove(end + 1) = i;
                break;
            end
        end
    end
    if ~isempty(remove)
        p(:, remove) = [];
        sigma(remove) = [];
        v(:, remove) = [];
        pv(remove) = [];
        ngroups = size(p, 2);
    end
    
    % break when v converges
    if all(size(v) == size(oldv)) && min(diag(oldv' * v)) > 0.999
        break;
    end    
    oldv = v;
    
end

% convert vanishing directions to [x y 1] form with (0,0) in upper-left of
% image
v = v';
v(:, 3) = v(:, 3) + 1E-4;
v = v ./ repmat(v(:, 3), 1, 3);
v(:, 1) = v(:, 1) + imsize(2) / imsize(1) / 2;
v(:, 2) = v(:, 2) + 1 / 2;

vpdata.v = v;
vpdata.vars = sigma;
vpdata.p = p;

end