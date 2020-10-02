%try to do proofs of disconnectedness of double-lobe region

%This is a farkas certificate of infeasibility
%there exists a function v where v <= 0 on X0, v >= 0 on X1, and v
%decreases along trajectories. Trajectories will never cross the surface
%v(t,x) = 0, and this level sets separates X0 and X1.

SOLVE = 1;
DRAW = 1;

if SOLVE
    opt = set_path_options;
    opt.x = sdpvar(2,1);
    opt.Tmax = 2;

    opt.scale = 0;

        order = 2;
    d =2*order;
    T = 2; %maximum time

    FEAS = 0;
    SET = 0;
    if FEAS
        opt.X0 = [1.25; -1];
        opt.X1 = [1.5; 0.5];

    else
        if SET
            opt.X0.ineq = 0.05 - (opt.x(1) + 0.75)^2 + (opt.x(2) + 0.25)^2;
    %         opt.X1.ineq = 0.05 - (opt.x(1) - 1.5)^2 + (opt.x(2)-0.5)^2;
            opt.X1 = [1.5; 0.5];
        else
            opt.X0 = [-0.75; 0];
%             opt.X0 = [-0.75; 0.5];
            opt.X1 = [1.5; 0.5];
%             opt.X1 = [1; 0.5];
%             opt.X1 = [1; -0.5];
%             opt.X1 = [0.5; 0.75]; %needs order 4
        end
    end

    opt.box = 0;

    %constraint set
    f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
    X.ineq = f(opt.x);
    X = fill_constraint(X);

    
    opt.time_indep = 1;
    
    opt.X = X;
    
    out = set_path_infeas_box_indep(opt, order);
end

if DRAW && out.farkas

    syms t [1 1]
    syms x [2 1]
    fy = f(x);
    vy = out.vval([t; x]);

    
    xl = [-2, 2];
    yl = [-2, 2];
   

    epsilon = 1e-5;


    figure(3)
    clf
    hold on
    Nlevel = 40;
    level_range = 2;
    fcontour(vy, [xl, yl], 'Fill', 'on', 'levellist', linspace(-level_range, level_range, Nlevel), 'DisplayName', 'v(x) contour')
    scatter(opt.X0(1),opt.X0(2), 100, 'ok', 'DisplayName', 'X0')
    scatter(opt.X1(1),opt.X1(2), 100, '*k', 'DisplayName', 'X1')
    fimplicit(fy == 0, [xl, yl], 'k', 'DisplayName','X')
    title('Farkas Infeasibility Certificate Contour')
%         legend('location', 'northeast')
    h = colorbar;
    ylabel(h, 'saturated v(x)')


    
    hold off
end
