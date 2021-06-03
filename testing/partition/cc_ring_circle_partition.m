SOLVE =1;
PLOT = 0;
SAMPLE = 0;
EVAL = 0;

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 1;

opt.scale = 1;
opt.verbose = 1;

order = 3;
d = 2*order;


%% support set (in unit ball)

R_ring_outer = 1;
R_ring_inner = 0.7;
R_circ = 0.5;

x2 = sum(opt.x.^2);
f = (R_ring_outer^2-x2)*(x2-R_ring_inner^2)*(x2-R_circ^2);

% f1 = -(x2-R_ring_inner^2)*(R_circ^2-x2);
f1 = (x2-R_ring_inner^2)*(x2-R_circ^2);
f2 = (R_ring_outer^2-x2);


%% initial sets
X0 = [0.8; 0];
X1 = [0; 0];


% X0 = [0.8 0.8; 0.0, 0.1];

Nth = 10;
th = linspace(0, 2*pi, Nth);
circ = [cos(th); sin(th)];

POINTS0 = 1;
POINTS1 = 1;

if POINTS0
% X0 = 0.8*circ;
% X0 = 0.9*circ + [0.05; 0.05];
% X0 = [0.8 -0.8; 0.0, 0.1];
X0 = [0.0; 0.8];
else
    
    R0 = 0.9;
    C0 = [0;0];
    X0 = struct('ineq', [], 'eq', 0.9^2 - sum((opt.x-C0).^2));
end

if POINTS1
% X1 = [0 0.3; 0 0];
X1 = [0;0.0];
end

%small circle
% X = struct('ineq', f, 'eq', []);
%  X = struct('ineq', [f1; f2], 'eq', []);

X_in = struct('ineq', R_circ^2-x2, 'eq', []);
X_out = struct('ineq', [(R_ring_outer^2-x2);(x2-R_ring_inner^2)], 'eq',[]);

X = {X_in, X_out};
%f1 = (x2-R_ring_inner^2)*(x2-R_circ^2);
% f2 = (R_ring_outer^2-x2);
% X1_infeas.ineq  = 0.2^2 - sum((opt.x - [1.5; 0.5]).^2);

% f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
% X.ineq = f(opt.x);
% X = fill_constraint(X);

% X_left.ineq  = R^2 - sum((opt.x - X0).^2);
% X_right.ineq = R^2 - sum((opt.x - X1).^2);

% X0_circ.ineq = R0^2 - sum((opt.x - X0).^2);
% X1_circ.ineq = R0^2 - sum((opt.x - X1).^2);

% X = {X_left, X_right};

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;
% opt.X0 = X0_circ;
% opt.X1 = X1_circ;

if SOLVE
% IM = set_manager(opt);
% out = IM.check_connected(d);

order_range = [1, 2];

spacing = [2; 2; 2];
% spacing= [3; 3; 3];

% spacing = [10; 4; 4];
SM = set_manager_partition(opt, spacing);

out = SM.climb_connected(order_range);

% out.farkas
end

if SAMPLE
%     Np = 100;
    [test, X_out_func]=constraint_eval(X_out, opt.x, X0);
    [test, X_in_func]=constraint_eval(X_in, opt.x, X0);
    X_func = @(pt) X_out_func(pt) || X_in_func(pt);
    
    supp_func = @(t,x) support_event(t, x, X_func);

    s_opt = set_sample_options;
    % x0 = @() [0;0];

    s_opt.x0 = @() X0;

    s_opt.Tmax = opt.Tmax;
    s_opt.dt = 0.1;
    s_opt.X_func = supp_func;

    Np = 15;

    out_sim=set_walk(Np, s_opt);

end

% 
% %TODO: write plotting code for ellipses
if PLOT &&  out.status == conn_status.Disconnected
    bplot = ring_circ_plotter(opt, out);
    
    bplot.circle = struct('R_ring_outer',R_ring_outer, 'R_ring_inner', R_ring_inner, 'R_circ', R_circ);
    bplot.set_plot();

end
