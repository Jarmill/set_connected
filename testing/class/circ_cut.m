SOLVE = 1;
PLOT = 1;


opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 2;


opt.scale = 0;

order = 4;
d = 2*order;

%points
X0 = [-1; 0];
X1 = [1; 0];

R = 1; %circle radius


% c_cut = [1; 0];
c_cut = [0; 0];
Rin = 0.4;
Rout = 0.5;


sr = sum((opt.x - c_cut).^2);
% X.ineq  = [R^2 - sum(opt.x.^2);  (sr-Rin)*(sr-Rout)];
% X.ineq  = [R^2 - sum(opt.x.^2);  opt.x(1)^2-Rin];

Xr.ineq = [R^2 - sum(opt.x.^2);  (Rin-sr)];
Xl.ineq = [R^2 - sum(opt.x.^2);  -(Rout-sr)];

X = {Xr; Xl};


opt.X = X;
opt.X0 = X0;
opt.X1 = X1;


if SOLVE
% IM = set_manager(opt);
% out = IM.check_connected(d);
% spacing= [1; 2; 2];
spacing= [2; 2; 2];
SM = set_manager_partition(opt, spacing);
out = SM.climb_connected(4);
% out = SM.climb_connected(3);
% out.farkas
end

% %TODO: write plotting code for ellipses
% if PLOT &&  out.status == conn_status.Disconnected
%     figure(1)
%     clf
%     syms tv [1 1];
%     syms xv [2 1];
%     hold on 
%     color0 = [0.4940, 0.1840, 0.5560];
%     color1 = [0.4660, 0.6740, 0.1880];
%     
%     v0 = out.func.v0(xv);
%     v1 = out.func.v1(xv);
%     limits = [-1.5, 1.5, -0.6, 0.6];
%     fimplicit(v0 == 1, limits,'DisplayName','v(0, x) = 1', 'Color', color0)
%     fimplicit(v1 == 0, limits,'DisplayName','v(T, x) = 0', 'Color', color1)
%     viscircles(X0', R, 'color', 'k')
%     viscircles(X1', R, 'color', 'k')
%     viscircles(X0', R0, 'LineStyle', '--', 'color', color0)
%     viscircles(X1', R0, 'LineStyle', '--', 'color', color1)
%     legend('location', 'south')
% %     lobe_plot = lobe_plotter(opt, out);
% %     lobe_plot.contour_2d();
% %     lobe_plot.contour_3d();
% 
% figure(2)
% clf
% hold on
% v = out.func.v([tv; xv]);
% fimplicit3(v==0, [0, 2, limits], 'MeshDensity', 120);
% 
% 
% end