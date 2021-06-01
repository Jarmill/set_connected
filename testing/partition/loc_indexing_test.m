box = [-1, 1;
    0, 2];

time_range = [0, 5];

spacing = [5,4,2];

loc = cell(spacing);

ind = 11;
loc{ind} = 'hi';
out = cell(size(spacing));
[out{:}] = ind2sub(spacing, ind)