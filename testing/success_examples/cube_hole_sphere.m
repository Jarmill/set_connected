% box hollowed by a sphere
opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(3,1);

% opt.Tmax = 1; %order 1 certificate
opt.Tmax = sqrt(3); %order 3 certificate

x = opt.x;
k = [1;2;3];
sx = sum(k.*(x.^2));

X0 = [0;0;0];

X1 = [0.9; 0.2; 0.2];

R = 0.5;
e = 0.1;

g = (sx - (R+e)^2) * ((R-e)^2 - sx);

X = struct('ineq', g);
pf = polyval_func(g, x);

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;

IM = set_manager(opt);
out = IM.climb_connected([1, 8]);