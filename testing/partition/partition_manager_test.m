%% set up the options structure

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
% opt.Tmax = 10;
options.Tmax = 151;
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

opt.box = [-1.75,2;
    -1.25,1.65];

%% define the manager

spacing= [1; 2; 1];

% spacing = [4; 2; 3];

% spacing = [3;1;1];

SM = set_manager_partition(opt, spacing);

% [loc, limits] = SM.make_locations(opt, spacing);

out = SM.check_connected(4);

% prog_mgr = SM.make_program(4)
% out = SM.solve_program(prog_mgr)