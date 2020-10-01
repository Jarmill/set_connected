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
%     order = 4;
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
%             opt.X0 = [-0.75; 0];
            opt.X0 = [-0.75; 0.5];
%             opt.X1 = [1.5; 0.5];
%             opt.X1 = [1; 0.5];
%             opt.X1 = [1; -0.5];
            opt.X1 = [0.5; 0.75]; %needs order 4
        end
    end

    opt.box = 0;

    %constraint set
    f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
    X.ineq = f(opt.x);
    X = fill_constraint(X);

    opt.X = X;

%     out = set_path_infeas(opt, order);
    out = set_path_infeas_box(opt, order);
end

if DRAW && out.farkas
    figure(2)
    clf
    syms t [1 1]
    syms x [2 1]
    fy = f(x);
    vy = out.vval([t; x]);

    hold on
    xl = [-2, 2];
    yl = [-2, 2];
    fimplicit(fy == 0, [xl, yl], 'DisplayName','X')

    epsilon = 1e-5;

    scatter(opt.X0(1),opt.X0(2), 100, 'ok', 'DisplayName', 'X0')
    scatter(opt.X1(1),opt.X1(2), 100, '*k', 'DisplayName', 'X1')
    
    vy0 = subs(vy, t, 0);
    color0 = [0.4940, 0.1840, 0.5560];
    colorT = [0.4660, 0.6740, 0.1880];
   
    
    fimplicit(vy0 == -1, [xl, yl], 'DisplayName','v(0, x0) <= -1' ,'Color', color0)
    fimplicit(vy0 == -epsilon, [xl, yl], ':', 'Color', color0, 'DisplayName','v(0, x0) < 0')    
    

    vyT = subs(vy, t, opt.Tmax);
    fimplicit(vyT == 1, [xl, yl], 'DisplayName','v(T, x1) >= 1', 'Color', colorT)
%     fimplicit(vyT == epsilon, [xl, yl], 'DisplayName','v(T, x1) > 0')
    
    fimplicit(vyT == epsilon, [xl, yl], ':', 'Color', colorT, 'DisplayName','v(T, x1) > 0')

    legend('location', 'northwest')
    
    title('Farkas Infeasibility Certificate')
    
    hold off
end
