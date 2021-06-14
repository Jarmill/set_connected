function [out_sim] = set_walk(Npoints, options)
%SET_WALK.M Summary of this function goes here
%   Detailed explanation goes here

out_sim = cell(Npoints, 1);

n = length(options.x0());

function ucurr = u_func()
    ucurr = 2*rand(n, 1)-1; 
    if options.u_boundary
        ucurr = ucurr/max(abs(ucurr));
%     else
    end
end



for i = 1:Npoints
    xcurr = options.x0();
    out_sim{i} = set_walk_single(xcurr, options.X_func, @u_func, options.Tmax, options.dt);
    
    
end


end

