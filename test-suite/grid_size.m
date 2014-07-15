np = (1:30) * 24;
s_grid = zeros(1, length(np));
for n = 1:length(np)
    cols = floor(sqrt(np(n)));
    rows = floor(np(n)/cols);
    s_grid(n) = rows * cols;
end

plot(np, s_grid, '-x')