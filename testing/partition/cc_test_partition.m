SOLVE = 1;
SAMPLE = 1;
EVAL = 0;
PLOT = 1;

%this example works without partitioning

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 20;
% opt.Tmax = 2;
% opt.Tmax = 151;
opt.verbose = 1;

opt.scale = 1;

X_bump = 0.2*[0 0 0 -1 -1 -1;
              0 -1 -2 0 -1 -2];
% X_bump = 0.2*[0 0 0 -1 -1 -1 -2 -2 -2;
%               0 -1 -2 0 -1 -2 0 -1 -2];
%           
% X0_infeas = [-0.75; 0.5];
X0_infeas = [-0.75; 0.3]+ X_bump;
% X0_infeas = [1.3; 0.9]; %quite feasible

%points
X1_infeas = [1.5; 0.5] + X_bump;

%small circle
% X1_infeas.ineq  = 0.2^2 - sum((opt.x - [1.5; 0.5]).^2);

f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
X.ineq = f(opt.x);
X = fill_constraint(X);


opt.X = X;

    opt.X0 = X0_infeas;
    opt.X1 = X1_infeas;

    order_range = [1, 3];
    
    opt.box = [-1.75,2;
    -1.25,1.65];

if SOLVE
% IM = set_manager(opt);
% out = IM.check_connected(d);

% spacing= [1; 2; 1];

spacing = [10; 4; 4];
SM = set_manager_partition(opt, spacing);

out = SM.climb_connected(order_range);
end


if SAMPLE
    %     Np = 100;
        [test, X_func]=constraint_eval(X, opt.x, [0;0]);

        supp_func = @(t,x) support_event(t, x, X_func);

    s_opt = set_sample_options;
    % x0 = @() [0;0];

    s_opt.x0 = @() opt.X0(:,randi( size(opt.X0, 2)));

    
    s_opt.u_boundary = 0;
    s_opt.Tmax = 10;
%     s_opt.dt = 0.05;
    s_opt.dt = 0.1;
    s_opt.X_func = supp_func;
    s_opt.nonneg_func = out.func.nonneg;

    Np = 30;

    out_sim=set_walk(Np, s_opt);

end


if SAMPLE || EVAL
   
    
    out_sim = set_traj_eval(out, out_sim);
end



if PLOT && out.status == conn_status.Disconnected
    lobe_plot = lobe_plotter(opt, out, out_sim);
    lobe_plot.contour_2d();
    lobe_plot.v_plot();
    lobe_plot.nonneg_traj();
    lobe_plot.nonneg_zeta();
%     for i = 1:Np
% %     out_sim = set_walk(x0(), X_func, @() u_func(2), Tmax, dt);
%         plot(out_sim{i}.x(1, :), out_sim{i}.x(2, :), 'c', 'HandleVisibility', 'off');
%     end 
    
%     lobe_plot.contour_3d();
end