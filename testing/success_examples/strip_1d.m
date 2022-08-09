%the same as twointervals_orig

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(1,1);
opt.Tmax = 1;

x = opt.x;

e = 0.1;
X0 = -(1+e)/2;
X1 = (1+e)/2;

% xleft  = 0.4;
% xright = 0.6;

g = -(1+x)*(-e-x)*(-e+x)*(1-x);

X = struct('ineq', g);

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;

IM = set_manager(opt);
out = IM.climb_connected([1, 8]);