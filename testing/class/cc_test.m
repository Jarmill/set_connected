opt = set_path_options;
opt.x = sdpvar(2,1);
opt.t = sdpvar(1, 1);
opt.Tmax = 2;


opt.scale = 0;

order = 4;
d = 2*order;

% X0_feas = [1.25; -1];
% X1_feas = [1.5; 0.5];

X0_infeas = [-0.75; 0.5];
X1_infeas = [1.5; 0.5];

f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
X.ineq = f(opt.x);
X = fill_constraint(X);

opt.X = X;
opt.X0 = X0_infeas;
opt.X1 = X1_infeas;
IM = infeas_manager(opt);

IM = IM.make_program(d);