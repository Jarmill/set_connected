%points in two different lobes of an elliptic curve (in 3d)

% ec = @(x) x.*(x+1).*(x-0.8).*(x+0.5).*(x-0.5);
ec = @(x) x(1).^3 -0.8*x(1).^2 +0.05;
ee = @(x,y,z) -((2*y.^2+0.5*z.^2) - ec(x));

fimplicit3(ee, [-1,1,-1,1,-1,1])

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(3,1);

% opt.Tmax = 1; %order 1 certificate
opt.Tmax = sqrt(3); %order 3 certificate

x = opt.x;

X0 = [0;0;0];

X1 = [0.9; 0.2; 0.2];

g = ee(opt.x(1), opt.x(2), opt.x(3));

X = struct('ineq', g);
pf = polyval_func(g, x);

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;

IM = set_manager(opt);
out = IM.climb_connected([1, 8]);