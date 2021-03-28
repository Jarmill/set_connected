opt = set_path_options;
opt.x = sdpvar(2,1);
opt.t = sdpvar(1, 1);
opt.Tmax = 2;


opt.scale = 0;

order = 3;
d = 2*order;

% X0_feas = [1.25; -1];
% X1_feas = [1.5; 0.5];

X1_bump =full(sparse([1 2], [1 3], [0.05, -0.05]));

X0_infeas = [-0.75; 0.5];
X1_infeas = [1.5; 0.5] + X1_bump;

f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
X.ineq = f(opt.x);
X = fill_constraint(X);

opt.X = X;
opt.X0 = X0_infeas;
opt.X1 = X1_infeas;
IM = infeas_manager(opt);

[cc_all, poly_var] = IM.make_program(d);


%TODO: 
%write solving code in class
opts = sdpsettings('solver', opt.solver);
opts.sos.model = 2;
[sol, monom, Gram, residual] = solvesos(cc_all.con, 0, opts, cc_all.coef);


%recovery code in class
%evaluation of functions
%visualization
[cv,mv] = coefficients(poly_var.v,[poly_var.t; poly_var.x]);
v_eval = value(cv)'*mv;

[cz, mz] = coefficients(poly_var.zeta,[poly_var.t; poly_var.x]);
zeta_eval = value(cz)*mz;