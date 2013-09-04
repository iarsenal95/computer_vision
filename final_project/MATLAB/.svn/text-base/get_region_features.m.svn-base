function rfeatures = get_region_features(image, imsegs, curmap, regind, spdata, vpdata)

% features are as follows:
% 1-12: mean abs oriented filter responses (ofr)
% 13: edginess - mean of all ofrs
% 14: index of largest ofr
% 15: dominance of largest response - max - median of ofr
% textons (would be 16-29) not often used 
%       16-27: mean abs texton filter response (tfr)
%       28: index of largest tfr
%       29: dominance of largest tfr (max - med)
% 16-18: red, green, blue means
% 19-21: hue, saturation, value means (rgb2hsv)
% 22-23: y location, x location means
% 24-29: hue - histogram (5 bins), entropy
% 30-33: sat - histogram (3 bins), entropy
% 34-37: location - 10% and 90% (percentiles) y and x 
% 38-39: edginess center y and x
% 40: number of super-pixels in segmentation
% 41: % of image area region takes up
% 42: number of sides in convex hull polygon
% 43: (num pixels) / (area of convex hull polygon)
% 44: whether the region is contiguous
% 45-60: vanishing point features (see vp2regionfeatures)
% 61-62: 10% and 90% y wrt horizon
% 63   : 0-> below horizon, 1-> straddles horizon, 2-> above horizon
% 64-82: older vanishing point features (see lines_to_vp_info3)
% discrete: 14 44 60 63

nclusters = length(regind);

rfeatures = zeros(nclusters, 82);

nangles = size(spdata.orientation, 2);

spfeatures = spdata.features; % superpixel features
pixelp = imsegs.npixels ./ sum(imsegs.npixels);
sporient = spdata.orientation;

height = size(image, 1);
width = size(image, 2);

hue = spdata.hue;
sat=  spdata.sat;
yim = 1 - repmat(((0 : height - 1) / (height - 1))', 1, width);
xim = repmat(((0 : width - 1) / (width - 1)), height, 1);

%% convert indices of superpixels to the indices of regions
nr = max(curmap);
rinds = cell(nr, 1);
for r = 1 : nr
    count = 0;
    rs = find(curmap == r);
    for k = 1 : length(rs)
        count = count + length(imsegs.seginds{rs(k)});
    end         
    rinds{r} = zeros(count, 1);
    
    count = 0;
    for k = 1 : length(rs) 
        rinds{r}(count + 1 : count + length(imsegs.seginds{rs(k)})) = imsegs.seginds{rs(k)};
        count = count + length(imsegs.seginds{rs(k)});
    end
end

%%
edgeim = spdata.edginess;

for c = 1 : nclusters
    
    spind = find(curmap == regind(c));
    cpixelp = pixelp(spind);
    
    % features consisting of mean block responses
    % mean orientation filter responses
    rorient = zeros(nangles, 1);
    for a = 1 : nangles
        rfeatures(c, a) = sum(spfeatures(spind, a) .* cpixelp);
        rorient(a) = sum(sporient(spind, a) .* cpixelp);        
    end
    % mean edginess (overall response)
    rfeatures(c, nangles + 1) = sum(spfeatures(spind, nangles + 1) .* cpixelp);
    % most dominant filter response
    [~, rfeatures(c, nangles + 2)] = max(rorient); 
    % dominance of largest filter (max - median)    
    rfeatures(c, nangles + 3) = max(rorient) - median(rorient);   

    nf = nangles + 3;
        
    % rgb means
    for b = 1 : 3
        rfeatures(c, nf + b) = sum(spfeatures(spind, nf + 2 + b) .* cpixelp);
    end
    
    % hsv
    rfeatures(c, nf + (4 : 6)) = rgb2hsv(rfeatures(c, nf + (1 : 3)));
    
    % y and x locs
    rfeatures(c, nf + 7) = sum(spfeatures(spind, nf + 7) .* cpixelp);
    rfeatures(c, nf + 8) = sum(spfeatures(spind, nf + 8) .* cpixelp);
    
    % get the values that pertain to this region
    rhue = hue(rinds{c});
    rsat = sat(rinds{c});
    rx = xim(rinds{c});
    ry = yim(rinds{c});
    redge = edgeim(rinds{c});
    npix = length(rinds{c});
    
    % hue and sat histograms
    nf = nf + 8;
    hue_histogram = (hist(rhue, 0.1 : 0.2 : 0.9) + 0.01) / (length(rhue) + 0.05);    
    rfeatures(c, nf + (1 : 5)) = hue_histogram / sum(hue_histogram);
    rfeatures(c, nf + 6) = -sum(hue_histogram .* log(hue_histogram)) / log(length(hue_histogram));
    sat_histogram = (hist(rsat, 0.167 : 0.333 : 0.833) + 0.01) / (length(rsat) + 0.03); 
    rfeatures(c, nf + (7 : 9)) = sat_histogram / sum(sat_histogram);
    rfeatures(c, nf + 10) = -sum(sat_histogram .* log(sat_histogram)) / log(length(sat_histogram));
    
    % location - 10% and 90% percentiles of y and x
    sorty = sort(ry);
    rfeatures(c, nf + 11) = sorty(ceil(npix / 10));
    rfeatures(c, nf + 12) = sorty(ceil(9 * npix / 10));
    sortx = sort(rx);
    rfeatures(c, nf + 13) = sortx(ceil(npix / 10));
    rfeatures(c, nf + 14) = sortx(ceil(9 * npix / 10));    
    
    % center of edginess y and x
    if rfeatures(c, nf + 11) == rfeatures(c, nf + 12)
        rfeatures(c, nf + 15) = 0.5;
    else
        center_y = sum(redge .* ry) / sum(redge);
        rfeatures(c, nf + 15) = sum((ry < center_y) .* ry) / sum(ry);
    end
    if rfeatures(c, nf + 13) == rfeatures(c, nf + 14)
        rfeatures(c, nf + 16) = 0.5;
    else
        center_x = sum(redge .* rx) / sum(redge);
        rfeatures(c, nf + 16) = sum((rx < center_x) .* rx) / sum(rx);
    end    

    % num superpixels, % of image area
    rfeatures(c, nf + 17) = length(spind);
    rfeatures(c, nf + 18) = npix / width / height;
    
    % polygon: num sides, area ratio
    try
        [polyi, polya] = convhull(rx, ry);
        rfeatures(c, nf + 19) = length(polyi) - 1;
        rfeatures(c, nf + 20) = npix / (polya * width * height);
    catch
        rfeatures(c, nf + (19 : 20)) = [4 0.75];
    end

    % whether contiguous
    rfeatures(c, nf + 21) = 1;
    for s = 1 : length(spind)
        isadj = imsegs.adjmat(setdiff(spind, spind(s)), spind(s));
        if sum(isadj) == 0
            rfeatures(c, nf + 21) = 0;         
            break;
        end
    end

    % vanishing point features
    nf = nf + 21;
    region_center = [sorty(ceil(npix / 2)) sortx(ceil(npix / 2))];  

    rbounds = [sortx(1) sortx(end) 1 - sorty(end) 1 - sorty(1)];
    rfeatures(c, nf + (1 : 16)) = vp2rfeatures(spind, vpdata, region_center, rbounds, imsegs);

    % y-location with respect to estimated horizon
    if ~isnan(vpdata.hpos)
        % location - 10% and 90% percentiles of y and x
        rfeatures(c, nf + 17) = rfeatures(c, 48) - (1 - vpdata.hpos); % bottom 10 pct wrt horizon
        rfeatures(c, nf + 18) = rfeatures(c, 49) - (1 - vpdata.hpos); % top 10 pct wrt horizon
        % 1 -> completely under horizon, 2-> straddles horizon, 3-> completely above horizon
        rfeatures(c, nf + 19) = (rfeatures(c, nf + 20) > 0) + (rfeatures(c, nf + 21) > 0) + 1;
    else % horizon was not estimated with high confidence
        rfeatures(c, nf + (17 : 18)) = rfeatures(c, 48 : 49) - 0.5;
        rfeatures(c, nf + 19) = 4;  % signifies no data-estimated horizon
    end
    
    rfeatures(c, nf + (20 : 38)) = lines2vpfeature(vpdata.spinfo(spind), vpdata.lines);  

end