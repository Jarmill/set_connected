SOLVE = 1;
PLOT = 1;

FEAS = 1;

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 2;


opt.scale = 0;


%small circle
% X1_infeas.ineq  = 0.2^2 - sum((opt.x - [1.5; 0.5]).^2);

% f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
a = 1;
% f = @(x) -(x(1).^2-a^2)*(x(1)-a).^2 + (x(2).^2-a^2).^2;
f = @(x, y) -(x.^2-a^2)*(x-a).^2 + (y.^2-a^2).^2;
%% 

X.ineq = f(opt.x(1), opt.x(2));
X = fill_constraint(X);

X0_feas = [0; 1];
X1_feas = [0; -1];

opt.X = X;
if FEAS
    opt.X0 = X0_feas;
    opt.X1 = X1_feas;
else
    opt.X0 = X0_infeas;
    opt.X1 = X1_infeas;
end
% opt.verbose = 1;
order_range = [1, 4];

if SOLVE
IM = set_manager(opt);
% out = IM.check_connected(d);
out = IM.climb_connected(order_range);
end

if PLOT
    bplot = bicuspid_plotter(opt, out.feas);
    bplot.set_plot();
end

% if PLOT && out.status == conn_status.Disconnected
%     lobe_plot = lobe_plotter(opt, out.infeas);
%     lobe_plot.contour_2d();
% %     lobe_plot.contour_3d();
% end