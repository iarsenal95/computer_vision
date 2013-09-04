function p = get_kde_likelihood(f, x, y)
n = length(x);
wx = x(2) - x(1);
indices = round((y - x(1)) / wx) + 1;
indices = min(max(indices, 1), n);
p = f(indices);
end