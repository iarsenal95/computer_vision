function hpos = estimate_horizon(vpdata, imsize)

v = vpdata.v;
vars = vpdata.vars;
p = vpdata.p;

% find vanishing points that occur within image height (or close)
ind = find(v(:, 2) > -0.25 & v(:, 2) < 1.25);

% do not consider vp with variance greater than varthresh
varthresh = 0.005;
ind(find(vars(ind) > varthresh)) = [];

% do not consider vp with fewer than 6 lines (expected)
memberthresh = 6;
nmembers = sum(p(:, ind), 1);
ind(find(nmembers < memberthresh)) = [];

v = v(ind, :);
vars = vars(ind);

% compute horizon position (if possible) and tilt (if possible)
if isempty(ind)
    hpos = NaN;
elseif ~isempty(ind)
    % change variables so variance is along y axis only
    v(:, 1) = v(:, 1) - 0.5 * imsize(2) / imsize(1);
    v(:, 2) = v(:, 2) - 0.5;
    for i = 1 : length(ind)
        vars(i) = vars(i) * v(i, 2) ^ 2 / (v(i, 1) ^ 2 + v(i, 2) ^ 2);    
    end    
    hpos = sum(v(:, 2) ./ vars') / sum(1 ./ vars');
    hpos = hpos + 0.5;
end

% plot the horizon
if exist('DISPLAY') && DISPLAY
    plot([1 imsize(2)], [hpos * imsize(1) hpos * imsize(1)], 'g', 'LineWidth', 2);
    drawnow;
end

end