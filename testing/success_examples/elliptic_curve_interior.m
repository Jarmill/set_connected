%points in two different lobes of an elliptic curve


opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 1;

x = opt.x;

R = 0.5;
e = 0.1;

% X0 = [0;0];
% X1 = [-1; 0.8];

% X0 = [-0.5; 0];
X0 = [-0.4; 0.1];
X1 = [0.8; 0.4];

% g = -(sum((M*(x-Xc)).^2) - (R+e)^2) * ((R-e)^2 - sum((M*(x-Xc)).^2));
g = -(x(2).^2 - x(1).^3 -0.8*x(1).^2 +0.05);

X = struct('ineq', g);
pf = polyval_func(g, x);

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;

IM = set_manager(opt);
out = IM.climb_connected([1, 8]);