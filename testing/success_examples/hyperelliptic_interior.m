% hc = @(x) x.*(x+1).*(x-3).*(x+2).*(x-2);
hc = @(x) x.*(x+1).*(x-0.8).*(x+0.5).*(x-0.5);
he = @(x,y) -(y.^2 - hc(x));
% fcontour(@(x,y) he(x,y), [-1,1,-1,1], 'levellist', 0)
% 
% figure(2)
% fsurf(@(x,y) he(x,y), [-1,1,-1,1])


%points in two different lobes of an elliptic curve


opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);

%order 1 at opt.Tmax = 1
opt.Tmax = 1;
%anything higher: indeterminate
% opt.Tmax = sqrt(2);

x = opt.x;

g = he(opt.x(1), opt.x(2));

X = struct('ineq', g);
pf = polyval_func(g, x);


X0 = [-0.85; -0.2];

X1 = [0.2; 0.1];

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;

IM = set_manager(opt);
out = IM.climb_connected([1, 8]);