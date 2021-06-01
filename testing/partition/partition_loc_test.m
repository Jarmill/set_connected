%% set up the options structure

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 2;
% opt.Tmax = 51;
opt.verbose = 1;

opt.scale = 1;

X_bump = 0.2*[0 0 0 -1 -1 -1;
              0 -1 -2 0 -1 -2];
X0_infeas = [-0.75; 0.3]+ X_bump;
% X0_infeas = [1.3; 0.9]; %quite feasible

%points
X1_infeas = [1.5; 0.5] + X_bump;

f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
X.ineq = f(opt.x);
X = fill_constraint(X);

opt.X = X;
opt.X0 = X0_infeas;
opt.X1 = X1_infeas;


%% define the location
time_range = [0, 2];
box = box_process(2, 1.5);
loc1 = set_location(opt, [0,1], box, 1);
loc2 = set_location(opt, [1,2], box, 2);

d = 4;
[poly_out1, coeff_out] = loc1.make_poly(d);
loc1.poly = poly_out1;

[poly_out2, coeff_out] = loc2.make_poly(d);
loc2.poly = poly_out2;

% pv = loc.get_adjacency_poly(0, 0);

loc1.next{1} = loc2;
loc2.prev{1} = loc1;

% con_adj = loc1.make_adjacency_con();

% [cons_loc2, coeff_loc2, nonneg, poly_out] = loc2.make_cons(d, poly_out2);

prog1 = loc1.make_program(d);
prog2 = loc2.make_program(d);