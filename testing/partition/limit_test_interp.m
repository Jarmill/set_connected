limits = {[             0 1 2 3 4 5];
    [-1.7500 0.1250 2];
    [  -1.2500 1.6500]};

spacing = cellfun(@length, limits);
dlim = cellfun(@(l) l(2) - l(1), limits);
start = cellfun(@(l) l(1), limits);

x = [ 1.33; 1; -1];

grid_handle = @(x) grid_ind(x, dlim, start, spacing);


ind = grid_ind(x, dlim, start, spacing);

% (x-start)./dlim;
% 
% sub = ceil((x-start)./dlim)';
% 
% sub(sub == 0) = 1;
% 
% sub = num2cell(sub);
% 
% ind = sub2ind(spacing, sub{:});

loc = cell(spacing');

loc(ind)