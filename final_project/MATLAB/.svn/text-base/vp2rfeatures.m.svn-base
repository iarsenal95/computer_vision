function features = vp2rfeatures(sinds, vpdata, rcenter, rbounds, imsegs)
% features
% 1 - 8 : line orientation histogram (normalized, weighted by length)
% 9     : entropy of line orientation histogram
% 10    : (num line pixels) / sqrt(area)
% 11    : (num line pixels with vertical vp membership) / sqrt(area)
% 12    : (num line pixels with horizontal vp membership) / sqrt(area)
% 13    : % of total line pixels with vertical membership
% 14    : x-pos of vp on horizon - center (0 if none)
% 15    : (y-pos of highest/lowest vertical vp - center) (0 if none)
% 16    : center wrt horizon vp {left, right, surrounding, far left/right, no vp on horizon}
% discrete features: 16
%

features = zeros(1, 16);

% maximum variance for a vp to be considered reliable
maxvpvar = 0.001;

% min distance to be considered "far" 
minfardist = 2.5; 

% max y-value for vp to be considered horizon candidate
maxhorzpos = 1.25;

% select relevant parts of data to this region
spinfo = vpdata.spinfo;
lines = vpdata.lines;
spinfo = spinfo(sinds);
lineinds = union([spinfo(:).lines], []);
lines = lines(lineinds, :);

% arrange vp info
nvp = size(vpdata.v, 1);
isvalidvp = ones(size(vpdata.vars));
vpdata.p = vpdata.p(lineinds, :);
[~, linevp] = max(vpdata.p, [], 2);

for i = 1 : length(isvalidvp)
    if isvalidvp(i) && sum(vpdata.p(:, i)) <= 3
        isvalidvp(i) = 0;
    end
end

% make vp relative to region center
vpdata.v = vpdata.v - repmat([rcenter 1], nvp, 1);

% more line and area statistics
nlines = size(lines, 1);
nlinepixels = sum([spinfo(:).nedges]);
linelength = sqrt((lines(:, 1) - lines(:, 2)) .^ 2 + (lines(:, 3) - lines(:, 4)) .^ 2);
rarea = sum(imsegs.npixels(sinds)) / prod(imsegs.imsize(1 : 2));

% orientation histogram and entropy
nbins = 8;
bw = pi / nbins;
binedges = bw / 2 : bw : pi - bw / 2;
lineangles = mod(lines(:, 5), pi); % range from 0 to pi
bins = zeros(nlines, 1);
for k = 1 : nlines
    bins(k) = sum(lineangles(k) > binedges);
end
bins = bins + nbins * (bins == 0);

linehist = zeros(nbins, 1);
for k = 1 : nbins
    linehist(k) = sum(linelength(find(bins == k)));
end
if sum(linehist) == 0
    features(1 : 8) = 0;
    features(9) = 1;
else
    linehist = linehist + 1;
    linehist = linehist / sum(linehist);
    features(1 : 8) = linehist;
    features(9) = -1 * sum(linehist .* log(linehist)) / log(nbins);
end

% 10    : (num line pixels) / sqrt(area)
features(10) = sum(nlinepixels) / sqrt(rarea);

% for vertical check that vanishing point has at least 75% confidence of
% being greater than 75 degrees
vertvp = find((abs(vpdata.v(:, 2)) > minfardist) & (vpdata.vars < maxvpvar)' & isvalidvp'); 
horzvp = find((abs(vpdata.v(:, 2)) < maxhorzpos) & isvalidvp' & (vpdata.vars < maxvpvar)'); 

% 11    : (num line pixels with vertical vp membership) / sqrt(area)
features(11) = 0;
for i = 1 : length(vertvp)
    vertlineinds = find(linevp == vertvp(i));
    features(11) = features(11) + sum(linelength(vertlineinds));
end
features(11) = features(11) / sqrt(rarea);

% 12    : (num line pixels with horizontal vp membership) / sqrt(area)
for i = 1 : length(horzvp)
    horzlineinds = find(linevp == horzvp(i));
    features(12) = features(12) + sum(linelength(horzlineinds));
end
features(12) = features(12) / sqrt(rarea);

% 13    : % of total line pixels with vertical membership
if linelength > 0
    features(13) = features(11) / ( sum(linelength) / sqrt(rarea));
end

% expected number of lines belonging to each vanishing point
vpcount = sum(vpdata.p, 1);

% 14    : x-pos of vp on horizon - center (0 if none)
[~, maxind] = max(vpcount(horzvp));
besthorzdist = vpdata.v(horzvp(maxind), 1);
if ~isempty(besthorzdist)
    features(14) = besthorzdist;
end % else features(14) = 0

% 15    : abs(y-pos of highest/lowest vertical vp - center) (0 if none)
[tmp, maxind] = max(vpcount(vertvp));
bestvertdist = vpdata.v(vertvp(maxind), 2);
if ~isempty(bestvertdist)
    features(15) = bestvertdist;
end % else features(15) = 0

% 16    : region pos wrt horizon vp {left, right, surrounding, far left/right, no vp on horizon}
tmppos = besthorzdist + rcenter(1); % return to image coordinates
if (isempty(besthorzdist)) % no vp on horizon
    features(16) = 5; 
elseif abs(besthorzdist) > minfardist % region far from vp (implies facing center)
    features(16) = 4;
elseif rbounds(2) < tmppos % region left of vp (implies facing right)
    features(16) = 1;
elseif rbounds(1) > tmppos % region right of vp (implies facing left)
    features(16) = 2; 
else  % region surrounds vp (implies mixed)
    features(16) = 3;
end

end