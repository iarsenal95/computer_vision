function gt_label
clear all;
close all;
clc;

%% add path
addpath(genpath(pwd));

%% load the input image
imname = 'images/vehicle10.jpg'; % change the imnames
global segname;
segname = [imname(1 : end - 4) '_sp.ppm'];

im = imread(imname);
imname = imname(1 : end - 3);
imname = [imname 'ppm'];
imwrite(im, imname);

%% segment the image
segcmd = [pwd '/segment/segment 0.8 100 100'];
system([segcmd ' ' imname ' ' segname]);

%% pre-process the segmented image
disp('processing...');
imsegs = process_segments(segname);
imsegs.imname = imname;
% remove the seginds field
%imsegs = rmfield(imsegs, 'seginds');
% add gvs_names field
imsegs.gvs_names = {'000', '090', 'sky'};

%% let user label the ground truth
global gvs;
gvs = zeros(imsegs.nseg, 1);

oriim = imread(imname);
imshow(oriim);

global segim;
segim = imread(segname);

% alpha blending
alpha = 0.2;
segim(:, :, 1) = segim(:, :, 1) * alpha + oriim(:, :, 1) * (1 - alpha);
segim(:, :, 2) = segim(:, :, 2) * alpha + oriim(:, :, 2) * (1 - alpha);
segim(:, :, 3) = segim(:, :, 3) * alpha + oriim(:, :, 3) * (1 - alpha);

global button_down;
button_down = 0;

% 1 - ground, 2 - vertical, 3 - sky
global label;
label = 0;
label_image(imsegs);

end

function label_image(imsegs)

global segim;
global label;
label = label + 1;

h = figure;
imshow(segim);
set(h, 'WindowButtonDownFcn', {@buttondown, imsegs});
set(h, 'WindowButtonUpFcn', {@buttonup});
set(h, 'WindowButtonMotionFcn', {@buttonmove, imsegs});
set(h, 'WindowKeyPressFcn', {@keypress, imsegs});
set(h, 'Pointer', 'fullcross');

switch label
    case 1
        label_name = 'Ground';
    case 2
        label_name = 'Vertical';
    case 3
        label_name = 'Sky';
end
fprintf('Now to label %s.\n', label_name);

end

function buttondown(~, ~, imsegs)
pt = get(gca, 'CurrentPoint');
x = uint16(pt(1, 1));
y = uint16(pt(1, 2));

width = imsegs.imsize(2);
height = imsegs.imsize(1);

global label;
global tlabel;
global segim;
global gvs;

% if right click, change the current label
type = get(gcbf, 'SelectionType');
if strcmp(type, 'alt')
    tlabel = mod(tlabel, 3) + 1;
else
    tlabel = label;
end

switch tlabel
    case 1
        label_name = 'Ground';
    case 2
        label_name = 'Vertical';
    case 3
        label_name = 'Sky';
end

if x > 0 && x < width && y > 0 && y < height
    % label the selected segment
    s = imsegs.segimage(y, x);
    fprintf('Label %d to %s\n', s, label_name);
    gvs(s) = tlabel;
    
    [a b] = find(imsegs.segimage == s);
    for i = 1 : length(a)
        segim(a(i), b(i), 1) = 0;
        segim(a(i), b(i), 2) = 0;
        segim(a(i), b(i), 3) = 0;
        segim(a(i), b(i), tlabel) = 255;
    end
    imshow(segim);
end

global button_down;
button_down = 1;

end

function buttonup(~, ~)
global button_down;
button_down = 0;
end
   
function buttonmove(~, ~, imsegs)
pt = get(gca, 'CurrentPoint');
x = uint16(pt(1, 1));
y = uint16(pt(1, 2));

width = imsegs.imsize(2);
height = imsegs.imsize(1);

global segim;
high_segim = segim;

global tlabel;
global gvs;
global button_down;

if x > 0 && x < width && y > 0 && y < height    
    % label the selected segment
    s = imsegs.segimage(y, x);
    
    if ~button_down
        [a b] = find(imsegs.segimage == s);
        for i = 1 : length(a)
            % highlight
            high_segim(a(i), b(i), 1) = min(high_segim(a(i), b(i), 1) * 1.5, 255);
            high_segim(a(i), b(i), 2) = min(high_segim(a(i), b(i), 2) * 1.5, 255);
            high_segim(a(i), b(i), 3) = min(high_segim(a(i), b(i), 3) * 1.5, 255);
        end
    else
        switch tlabel
            case 1
                label_name = 'Ground';
            case 2
                label_name = 'Vertical';
            case 3
                label_name = 'Sky';
        end
        
        fprintf('Label %d to %s\n', s, label_name);
        gvs(s) = tlabel;

        [a b] = find(imsegs.segimage == s);
        for i = 1 : length(a)
            segim(a(i), b(i), 1) = 0;
            segim(a(i), b(i), 2) = 0;
            segim(a(i), b(i), 3) = 0;
            segim(a(i), b(i), tlabel) = 255;
        end
    end
    imshow(high_segim);
end

end

function keypress(~, ~, imsegs)
char = get(gcf, 'CurrentCharacter');

global segname;
global label;
global gvs;
if char == 13
    close;
    if label < 3
        label_image(imsegs);
    else
        % delete the superpixel image
        delete(segname);
        
        % save the file
        % add gvs field
        imsegs.gvs = gvs;
        save([imsegs.imname '.mat'], 'imsegs');
        disp('Done!');
        close;
    end
end
end