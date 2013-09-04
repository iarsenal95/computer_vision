function [lines spdata] = find_long_edges(grayim, imsegs)

% reduce the noise and compute the gradient along x and y axis
[dx dy] = gradient(conv2(grayim, fspecial('gaussian', 7, 1.5), 'same'));

% find edges using canny detector
[cannyim ~] = edge(grayim, 'canny');

% remove border edges
cannyim([1 2 end - 1 end], :) = 0;
cannyim(:, [1 2 end - 1 end]) = 0;
[h w] = size(cannyim);

% figure;
% imshow(cannyim);

% get points on the edges
ind = find(cannyim > 0);

% get the gradient of these points
indx = dx(ind);
indy = dy(ind);
inda = atan(indy ./ (indx + 1E-10));

% inda ranges from 1 to ndir with bin centered around pi / 2
ndir = 8; % 8 directions
inda = ceil(mod(inda/ pi * ndir - 0.5, ndir));

% get the indices of edges in each direction
dirs = cell(1, ndir);
for i = 1 : ndir
    dirs{i} = uint32(ind(find(inda == i)));
end

lines = zeros(2000, 6);
used = zeros(size(cannyim));

nline = 0;

nspdata = length(imsegs.npixels);
spdata = repmat(struct('lines', zeros(5, 1), 'nedges', 0), nspdata, 1);
bcount = zeros(nspdata, 1);

if exist('DISPLAY') && DISPLAY
    figure;
    imshow(grayim); hold on;
end

% for each direction
for k = 1 : ndir
    
    nind = 0;
    % for 3 adjacent directions, pick the unused edge points
    for m = (k - 1) : k + 1
        nind = nind + sum(~used(dirs{mod(m - 1, ndir) + 1}));
    end
    
    ind = zeros(nind, 1);
    dirim = zeros(size(cannyim));
  
    count = 0;
    for m = (k - 1) : k + 1
        m2 = mod(m - 1, ndir) + 1;
        tind = dirs{m2}(find(~used(dirs{m2})));
        tmpcount = length(tind);
        ind(count + 1 : count + tmpcount) = tind;
        count = count + tmpcount;
    end
    dirim(ind) = 1;
        
    % compute the 8-connected neighbor regions
    [tmpL, nedge] = bwlabel(dirim, 8); 
    
    % get the number of pixels in each edge
    sedge = zeros(nedge, 1);
    edges = repmat(struct('ind', zeros(200, 1)), nedge, 1);
    for i = 1  :length(ind)
        id = tmpL(ind(i));
        sedge(id) = sedge(id) + 1;
        edges(id).ind(sedge(id)) = ind(i);
    end          
    for i = 1 : nedge
        edges(i).ind = edges(i).ind(1 : sedge(i));
    end        
       
    % get the endpoints of the long edges and an image of the long edges
    % using the method mentioned in "Video Compass"
    minlen =  sqrt(h ^ 2 + w ^ 2) * 0.02;
    for id = 1 : nedge
        if sedge(id) > minlen
            y = mod(edges(id).ind - 1, h) + 1;
            x = floor((edges(id).ind - 1) / h) + 1;
            
            mean_x = mean(x);
            mean_y = mean(y);
            zmx = (x - mean_x);
            zmy = (y - mean_y);
            D = [sum(zmx .^ 2) sum(zmx .* zmy); sum(zmx .* zmy) sum(zmy .^ 2)];
            [v, lambda] = eig(D);
            theta = atan2(v(2, 2) , v(1, 2));
            
            % the quality of the line is given by the ratio of the two
            % eigenvalues of D
            if lambda(1, 1) > 0
                conf = lambda(2, 2) / lambda(1, 1);
            else
                conf = 100000;
            end
            
            if conf >= 400 
                nline = nline + 1;
                
                % eliminate the used points
                used(edges(id).ind) = 1;
                
                bi = double(imsegs.segimage(edges(id).ind));
                [~, gn] = grp2idx(bi);
                for s = 1 : length(bi)
                    if bi(s) > 0
                        spdata(bi(s)).nedges = spdata(bi(s)).nedges + 1;
                    end
                end
                for s = 1:length(gn)
                    tmpbi = str2num(gn{s});
                    if tmpbi > 0
                        bcount(tmpbi) = bcount(tmpbi) + 1;
                        spdata(tmpbi).lines(bcount(tmpbi)) = nline;                
                    end
                end                
                
                r = sqrt((max(x) - min(x)) ^ 2 + (max(y) - min(y)) ^ 2);
                x1 = mean_x - cos(theta) * r / 2;
                x2 = mean_x + cos(theta) * r / 2;
                y1 = mean_y - sin(theta) * r / 2;
                y2 = mean_y + sin(theta) * r / 2;            
                
                r = mean_x * cos(theta) + mean_y * sin(theta);
                
                lines(nline, 1 : 6) = [x1 x2 y1 y2 theta r];
                
                if exist('DISPLAY') && DISPLAY
                    plot([x1 x2], [y1 y2], 'r');
                end
            end
        end
    end

end

for k = 1 : length(spdata)
    spdata(k).lines = spdata(k).lines(1 : bcount(k))';
end
lines = lines(1 : nline, :);

% then normalize the lines
if size(lines, 1) > 0
    lines(:, 1 : 2) = (lines(:, 1 : 2) - w / 2) / h;
    lines(:, 3 : 4) = (lines(:, 3 : 4) - h / 2) / h;
    lines(:, 6) = mean(lines(:, 1 : 2), 2) .* cos(lines(:, 5)) + mean(lines(:, 3 : 4), 2) .* sin(lines(:, 5));
end

end