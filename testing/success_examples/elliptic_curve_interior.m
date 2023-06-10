%points in two different lobes of an elliptic curve

SOLVE = 1;
PLOT = 1;

if SOLVE
opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);

opt.Tmax = 1; %order 1 certificate
opt.Tmax = sqrt(2); %order 3 certificate

x = opt.x;

X0 = [-0.4; 0.1];
X1 = [0.8; 0.4];

g = -(x(2).^2 - x(1).^3 -0.8*x(1).^2 +0.05);

X = struct('ineq', g);
pf = polyval_func(g, x);

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;

IM = set_manager(opt);
out = IM.climb_connected([1, 8]);

end

if PLOT

end