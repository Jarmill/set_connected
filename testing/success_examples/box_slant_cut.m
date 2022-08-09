%cutting the box with a slanted strip

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 1;

x = opt.x;

e = 0.1;
X0 = -((1+e)/2)*[1; 0];
X1 = ((1+e)/2)*[1; 0.5];

g = (e^2 - sum(x)^2);

X = struct('ineq', g);

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;

IM = set_manager(opt);
out = IM.climb_connected([1, 8]);