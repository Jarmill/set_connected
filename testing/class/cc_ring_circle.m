SOLVE = 1;
PLOT = 1;


opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 3;

opt.scale = 0;
opt.verbose = 1;

order = 5;
d = 2*order;


%% support set (in unit ball)

R_ring_outer = 1;
R_ring_inner = 0.7;
R_circ = 0.5;

x2 = sum(opt.x.^2);
% f = (R_ring_outer^2-x2)*(x2-R_ring_inner^2)*(R_circ^2-x2);

f1 = (x2-R_ring_inner^2)*(R_circ^2-x2);
f2 = (R_ring_outer^2-x2);


%% initial sets
X0 = [0.8; 0];
X1 = [0; 0];


% X0 = [0.8 0.8; 0.0, 0.1];

Nth = 10;
th = linspace(0, 2*pi, Nth);
circ = [cos(th); sin(th)];

POINTS0 = 0;
POINTS1 = 1;

if POINTS0
% X0 = 0.8*circ;
% X0 = 0.9*circ + [0.05; 0.05];
% X0 = [0.8 -0.8; 0.0, 0.1];
X0 = [0.8; 0.0];
else
    
    R0 = 0.9;
    C0 = [0;0];
    X0 = struct('ineq', [], 'eq', 0.9^2 - sum((opt.x-C0).^2));
end

if POINTS1
X1 = [0 0.3; 0 0];
end

%small circle
% X = struct('ineq', f, 'eq', []);
 X = struct('ineq', [f1; f2], 'eq', []);


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
IM = set_manager(opt);
out = IM.check_connected(d);
% out.farkas
end
% 
% %TODO: write plotting code for ellipses
if PLOT &&  out.status == conn_status.Disconnected
    figure(1)
    clf
    syms tv [1 1];
    syms xv [2 1];
    hold on 
    color0 = [0.4940, 0.1840, 0.5560];
    color1 = [0.4660, 0.6740, 0.1880];
%     

th_fine = linspace(0, 2*pi, 200);
    circ_fine = [cos(th_fine); sin(th_fine)];
    

    color_int = 0.8*[1,1,1];
    patch(R_ring_outer*circ_fine(1, :), R_ring_outer*circ_fine(2, :), color_int, 'LineWidth', 3, 'DisplayName', 'X')
    patch(R_ring_inner*circ_fine(1, :), R_ring_inner*circ_fine(2, :), 'w', 'LineWidth', 3, 'HandleVisibility', 'off')
    patch(R_circ*circ_fine(1, :), R_circ*circ_fine(2, :), color_int, 'LineWidth', 3, 'HandleVisibility', 'off')

    
    if POINTS0
        scatter(X0(1, :), X0(2, :), 100, 'ok', 'DisplayName', 'X0');
    else
        circ0 = R0*circ_fine + C0;
        plot(circ0(1, :), circ0(2, :), '--k', 'LineWidth', 3, 'DisplayName', 'X0');
    end
    scatter(X1(1, :), X1(2, :), 100, '*k', 'DisplayName', 'X1');

    
    v0 = out.func.v0(xv);
    v1 = out.func.v1(xv);
    limits = 1.5*[-1,1,-1,1];
    fimplicit(v0 == 1, limits,'DisplayName','v(0, x) = 1', 'Color', color0, 'LineWidth', 3);
    fimplicit(v1 == 0, limits,'DisplayName','v(T, x) = 0', 'Color', color1, 'LineWidth', 3);
    
    
    
    
    
    

    axis square;
    
%     viscircles(X1', R, 'color', 'k')
%     viscircles(X0', R0, 'LineStyle', '--', 'color', color0)
%     viscircles(X1', R0, 'LineStyle', '--', 'color', color1)
    legend('location', 'northeast', 'FontSize', 12)
    
    xlabel('x_1', 'FontSize', 12)
    ylabel('x_2', 'FontSize', 12)
    title(['Proof of Disconnectedness at order ', num2str(order)], 'FontSize', 16)
    
% %     lobe_plot = lobe_plotter(opt, out);
% %     lobe_plot.contour_2d();
% %     lobe_plot.contour_3d();
% 
% figure(2)
% clf
% hold on
% v = out.infeas.func.v([tv; xv]);
% fimplicit3(v==0, [0, 2, limits], 'MeshDensity', 120);
% 
% 
end
