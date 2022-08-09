%A circular cut out of the side

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
% opt.Tmax = 1;
opt.Tmax = 2;

x = opt.x;

R = 0.5;
e = 0.1;

% X0 = [0;0];
% X1 = [-1; 0.8];

X0 = [1; 0];
X1 = [-1; 1];
Xc = X0;

g = -(sum((x-Xc).^2) - (R+e)^2) * ((R-e)^2 - sum((x-Xc).^2));

X = struct('ineq', g);
pf = polyval_func(g, x);

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;

IM = set_manager(opt);
out = IM.climb_connected([1, 8]);