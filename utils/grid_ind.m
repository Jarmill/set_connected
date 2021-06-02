function [ind] = grid_ind(x, arg2, arg3, arg4)
%GRID_IND find the index of point x in a grid

%inputs:
%   x:     point in space
%2 arguments
%   arg2:  limits to the grid

%4 arguments
%   arg2:  grid spacing vector for each dimension
%   arg3:  start of each grid dimension
%   arg4:  number of elements in each direction
if nargin == 2
    limits = arg2;
    dlim = cellfun(@(l) l(2) - l(1), limits);
    start = cellfun(@(l) l(1), limits);
    spacing = size(limits);
else
    dlim = arg2;
    start = arg3;
    spacing = arg4;
end
    
%find the subscript of the cell in which x resides
sub = ceil((x-start)./dlim);

sub(sub == 0) = 1;


sub(sub >= spacing) = spacing(sub >= spacing);

sub = num2cell(sub)';
%output the index
ind = sub2ind(spacing, sub{:});
end

