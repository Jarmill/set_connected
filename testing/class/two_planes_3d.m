SOLVE = 0;
PLOT = 1;


opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(3,1);
opt.Tmax = 3;

opt.scale = 0;
opt.verbose = 1;

order = 2;
d = 2*order;


%% support set (in unit ball)

z_width = 0.5;


%% initial sets

X0 = [0;0; z_width];
X1 = [0;0; -z_width];
%small circle
% X = struct('ineq', f, 'eq', []);
 X = {struct('ineq', 1-opt.x(1:2).^2, 'eq', opt.x(3) - z_width),
     struct('ineq', 1-opt.x(1:2).^2, 'eq', opt.x(3) + z_width)};

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
    syms xv [3 1];
    hold on 
    color0 = [0.4940, 0.1840, 0.5560];
    color1 = [0.4660, 0.6740, 0.1880];
%     
    patch([1,-1,-1,1,1],[1,1,-1,-1,1],z_width*ones(5, 1), 'k', 'facealpha', 0.5)
    patch([1,-1,-1,1,1],[1,1,-1,-1,1],-z_width*ones(5, 1), 'k', 'facealpha', 0.5)
    color_int = 0.8*[1,1,1];
%     patch(R_ring_outer*circ_fine(1, :), R_ring_outer*circ_fine(2, :), color_int, 'LineWidth', 3, 'DisplayName', 'X')
%     patch(R_ring_inner*circ_fine(1, :), R_ring_inner*circ_fine(2, :), 'w', 'LineWidth', 3, 'HandleVisibility', 'off')
%     patch(R_circ*circ_fine(1, :), R_circ*circ_fine(2, :), color_int, 'LineWidth', 3, 'HandleVisibility', 'off')

    

    
%     if POINTS0
        scatter3(X0(1, :), X0(2, :), X0(3, :), 100, 'ok', 'DisplayName', 'X0');
%     else
%         circ0 = R0*circ_fine + C0;
%         plot(circ0(1, :), circ0(2, :), '--k', 'LineWidth', 3, 'DisplayName', 'X0');
%     end
    scatter3(X1(1, :), X1(2, :),X1(3, :), 100, '*k', 'DisplayName', 'X1');

    
    v0 = out.func.v0(xv);
    v1 = out.func.v1(xv);
    limits = [1.5*[-1,1,-1,1], 0.75*[-1,1]];
    fimplicit3(v0 == 1, limits,'DisplayName','v(0, x) = 1', 'FaceColor', color0, 'LineWidth', 3);
    fimplicit3(v1 == 0, limits,'DisplayName','v(T, x) = 0', 'FaceColor', color1, 'LineWidth', 3);
    
    
%     th_fine = linspace(0, 2*pi, 200);
%     circ_fine = [cos(th_fine); sin(th_fine)];
    
    
    
    

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