SOLVE =0;
PLOT = 1;
SAMPLE = 1;
EVAL = 1;

%ring-circle test exhibits circular symmetry, depends only on the radius
%solve a 1d set connectedness program based on the radius 

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(1,1);
opt.Tmax = 2;
opt.box = [0, 1];

opt.scale = 1;
opt.verbose = 1;



%% support set (in unit ball)

R_ring_outer = 1;
R_ring_inner = 0.7;
R_circ = 0.6;

x = opt.x;

% x2 = sum(opt.x.^2);
% f = (R_ring_outer^2-x2)*(x2-R_ring_inner^2)*(x2-R_circ^2);

% f1 = -(x2-R_ring_inner^2)*(R_circ^2-x2);
% f1 = (x2-R_ring_inner^2)*(x2-R_circ^2);
% f2 = (R_ring_outer^2-x2);


%% initial sets
X0 = [0.8; 0];
X1 = [0; 0];

POINTS0 = 1;
POINTS1 = 1;


X0 = [0.8];
X1 = [0.3];

X_in = struct('ineq', [R_circ - x; x] , 'eq', []);
X_out = struct('ineq', [R_ring_outer-x; x-R_ring_inner], 'eq',[]);

% X_in = struct('ineq', [(R_circ - x)* x] , 'eq', []);
% X_out = struct('ineq', [(R_ring_outer-x)*( x-R_ring_inner)], 'eq',[]);


X = {X_in, X_out};

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;
opt.epsilon = 0.01;

if SOLVE
% IM = set_manager(opt);

% spacing = [1; 1];

% spacing = [1; 2];
% spacing = [4;2];

spacing = [5; 6];
% spacing= [3; 3; 3];

% spacing = [10; 4; 4];
SM = set_manager_partition(opt, spacing);


order = 5;
% order = 4;
% order = 3;
d = 2*order;

out = SM.check_connected(d);
% out.farkas
end


if SAMPLE
%     Np = 100;


    xx = sdpvar(2, 1);
    x2= sum(xx.^2);
    f1 = (R_circ^2-x2);
    f2 = [(R_ring_outer^2-x2);(x2-R_ring_inner^2)];

    [test, X_out_func]=constraint_eval(struct('ineq', f1), xx, [0;0]);
    [test, X_in_func]=constraint_eval(struct('ineq', f2), xx, [0;0]);
    X_func = @(pt) X_out_func(pt) || X_in_func(pt);
    
    supp_func = @(t,x) support_event(t, x, X_func);

    s_opt = set_sample_options;
    % x0 = @() [0;0];

    s_opt.x0 = @() sphere_sample(1, 2)'*X0;

    s_opt.Tmax = opt.Tmax;
    s_opt.dt = 0.1;
    s_opt.X_func = supp_func;

%     Np = 200;
    Np = 300;

    out_sim=set_walk(Np, s_opt);

end

if (SAMPLE || EVAL) && (out.status == conn_status.Disconnected)
   
    out_sim_rad = out_sim;
    for i = 1:length(out_sim_rad)
        out_sim_rad{i}.x_traj = sqrt(sum(out_sim_rad{i}.x_traj.^2, 2));
        out_sim_rad{i}.x      = sqrt(sum(out_sim_rad{i}.x.^2, 1));
    end
    out_sim_rad = set_traj_eval(out, out_sim_rad);
    for i = 1:length(out_sim)
        out_sim{i}.v = out_sim_rad{i}.v;
        out_sim{i}.nonneg = out_sim_rad{i}.nonneg;
    end
end

% 
% %TODO: write plotting code for ellipses
if PLOT &&  out.status == conn_status.Disconnected
%     bplot = ring_circ_plotter(opt, out);
    bplot = ring_circ_plotter(opt, out, out_sim);

%     bplot.rad0 = R_circ;
%     bplot.opt.X0 = R_circ;
    
    bplot.circle = struct('R_ring_outer',R_ring_outer, 'R_ring_inner', R_ring_inner, 'R_circ', R_circ);
    bplot.set_plot();
    bplot.v_plot();
    bplot.nonneg_traj();
    bplot.nonneg_zeta();
        
    bplot.contour_2d();
end
