% only works when after calling main
im = imread(imname);
segim = imread(segname);

%% show some constellations
Nc = 5;
ind = find(npart == Nc);
ind = ind(1);
reds = zeros(Nc, 1);
greens = zeros(Nc, 1);
blues = zeros(Nc, 1);
for i = 1 : Nc
    reds(i) = unidrnd(256) - 1;
    greens(i) = unidrnd(256) - 1;
    blues(i) = unidrnd(256) - 1;
end

consim = im;
for i = 1 : size(segim, 1)
    for j = 1 : size(segim, 2)
        c = imsegs.segimage(i, j);
        c = maps(c, ind);
        consim(i, j, 1) = reds(c);
        consim(i, j, 2) = greens(c);
        consim(i, j, 3) = blues(c);
    end
end

imshow(consim);

%% show the final label
Nc = 3;
reds = zeros(Nc, 1);
greens = zeros(Nc, 1);
blues = zeros(Nc, 1);
for i = 1 : Nc
    reds(i) = unidrnd(256) - 1;
    greens(i) = unidrnd(256) - 1;
    blues(i) = unidrnd(256) - 1;
end

consim = im;
for i = 1 : size(segim, 1)
    for j = 1 : size(segim, 2)
        c = imsegs.segimage(i, j);
        label = labels.vert_labels(c);
        if strcmp(label, '000')
            c = 1;
        elseif strcmp(label, '090')
            c = 2;
        else
            c = 3;
        end
        consim(i, j, 1) = reds(c);
        consim(i, j, 2) = greens(c);
        consim(i, j, 3) = blues(c);
    end
end

figure;
imshow(consim);