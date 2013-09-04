function spmeans = get_spmeans(imseg, values)

spmeans = zeros(imseg.nseg, 1);
count = zeros(imseg.nseg, 1);
segimage = imseg.segimage;
nelements = numel(segimage);
for n = 1 : nelements
    s = segimage(n);
    if s~=0
        spmeans(s) = spmeans(s) + values(n);
        count(s) = count(s) + 1;
    end
end

spmeans = spmeans ./ (count + (count == 0));

end

    