%moon test 
%there was a bug in the code, X was set to {} before adding circles
%now the certificate fails. Will try partitioning


SOLVE = 0;
PLOT = 1;
SAMPLE = 1;
EVAL = 0;
FEAS = 0;

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 2;
opt.box = [-1,1;0,1];
opt.epsilon = 0.01;


opt.scale = 1;
opt.verbose = 1;

%moon set
inner_rad = 0.7;
inner_x = 0.35;

%use only upper half-plane

X_moon = struct;
X_moon.ineq = [1 - opt.x(1)^2 - opt.x(2)^2; 
                (opt.x(1)-inner_x)^2 + opt.x(2)^2 - inner_rad^2;
                opt.x(2)];

X = {X_moon};

%circles
MULTIPLE_CIRCLES = 0;

if MULTIPLE_CIRCLES
    circle_rad = [0.4; 0.3; 0.3];
    circle_center = [0.4, 0;
                 -1, 1;
                 -1, -1]';

%     circle_rad = [0.4; 0.3];
% circle_center = [0.4, 0;
%                  -1, 1]';

else    
    circle_rad = [0.4];
    circle_center = [0.4, 0]';
end

X_circ = {};
for i = 1:length(circle_rad)
%     X_circ = struct;
    rx = circle_rad(i);
    cx = circle_center(:, i);
    X_circ.ineq = [rx^2 - (opt.x(1)-cx(1))^2 - (opt.x(2)-cx(2))^2;
                    opt.x(2)];
    X = [X; {X_circ}];
end



X0_infeas = [0; 0.9];
X1_infeas = [0.4; 0];

opt.X = X;
opt.X0 = X0_infeas;
opt.X1 = X1_infeas;
order_range = [1, 4];

if SOLVE
% IM = set_manager(opt);
spacing = [4;4;4];

SM = set_manager_partition(opt, spacing);
order = 4;
d = 2*order;
out = SM.check_connected(d);
% out = IM.climb_connected(order_range);
end

if SAMPLE
    [test, X_moon_func]=constraint_eval(X_moon, opt.x, [0;0]);
    [test, X_circ_func]=constraint_eval(X_circ, opt.x, [0;0]);
    X_func = @(pt) X_moon_func(pt) || X_circ_func(pt);

    supp_func = @(t,x) support_event(t, x, X_func);

    s_opt = set_sample_options;
    % x0 = @() [0;0];

    s_opt.x0 = @() opt.X0();

    s_opt.Tmax = opt.Tmax;
    s_opt.dt = 0.04;
    s_opt.X_func = supp_func;

    Np = 300;

    out_sim=set_walk(Np, s_opt);
end


if (SAMPLE || EVAL) && (out.status == conn_status.Disconnected)
  
    out_sim = set_traj_eval(out, out_sim);

end

if PLOT
    bplot = moon_half_plotter(opt, out, out_sim);
    bplot.inner_x = inner_x;
    bplot.inner_rad = inner_rad;
    bplot.circle = struct('center', circle_center, 'rad', circle_rad);
    bplot.set_plot();
    bplot.traj_2d();
    
    bplot.v_plot();
    bplot.nonneg_traj();
    bplot.nonneg_zeta();
    bplot.contour_2d(4);
end