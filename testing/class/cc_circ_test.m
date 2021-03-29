SOLVE = 0;
PLOT = 1;


opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 2;


opt.scale = 0;

order = 3;
d = 2*order;

%points
X0 = [-1; 0];
X1 = [1; 0];

R = 0.5; %circle radius
R0 = 0.2;

%small circle
% X1_infeas.ineq  = 0.2^2 - sum((opt.x - [1.5; 0.5]).^2);

% f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
% X.ineq = f(opt.x);
% X = fill_constraint(X);

X_left.ineq  = R^2 - sum((opt.x - X0).^2);
X_right.ineq = R^2 - sum((opt.x - X1).^2);

X0_circ.ineq = R0^2 - sum((opt.x - X0).^2);
X1_circ.ineq = R0^2 - sum((opt.x - X1).^2);

X = {X_left, X_right};

opt.X = X;
% opt.X0 = X0;
% opt.X1 = X1;
opt.X0 = X0_circ;
opt.X1 = X1_circ;

if SOLVE
IM = set_infeas_manager(opt);
out = IM.infeas(d);
out.farkas
end

%TODO: write plotting code for ellipses
if PLOT && out.farkas
    figure(1)
    clf
    syms tv [1 1];
    syms xv [2 1];
    hold on 
    color0 = [0.4940, 0.1840, 0.5560];
    color1 = [0.4660, 0.6740, 0.1880];
    
    v0 = out.func.v0(xv);
    v1 = out.func.v1(xv);
    limits = [-1.5, 1.5, -0.6, 0.6];
    fimplicit(v0 == 1, limits,'DisplayName','v(0, x) = 1', 'Color', color0)
    fimplicit(v1 == -1, limits,'DisplayName','v(T, x) = -1', 'Color', color1)
    viscircles(X0', R, 'color', 'k')
    viscircles(X1', R, 'color', 'k')
    viscircles(X0', R0, 'color', 'k')
    viscircles(X1', R0, 'color', 'k')
%     lobe_plot = lobe_plotter(opt, out);
%     lobe_plot.contour_2d();
%     lobe_plot.contour_3d();
end