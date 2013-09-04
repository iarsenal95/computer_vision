function labels = get_region_labels(gt, curmap, curnsegs)
    labels = zeros(curnsegs, 1);
    for i = 1 : curnsegs
        ind = find(curmap == i);
        gtlabel = gt(ind);
        label = 4;
        for j = 1 : 3
            if length(find(gtlabel == j)) / nnz(gtlabel) > 0.9
                label = j;
                break;
            end
        end
        labels(i) = label;
    end
end