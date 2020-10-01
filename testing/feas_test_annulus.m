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

    opt.scale = 1;


    order = 2;
    d =2*order;
    T = 3; %maximum time

    FEAS = 1;
    SEP = 0;
    
    if FEAS
        opt.X0 = [1; -1];
        opt.X1 = [1; 1];
    else
        opt.X0 = [1; -1];
        opt.X1 = [-1; -1];
    end

    %constraint set
%     Rinner = 1.4;
    Rinner = 1;
    Router = 1.8;
    Absx = 0.5;

    x = opt.x;
    X = struct;
    if SEP
        X.ineq = [x'*x - Rinner^2; Router^2 - x'*x];
    else
        X.ineq = [x'*x - Rinner^2; Router^2 - x'*x; x(1)^2 - Absx^2];
    end
    X = fill_constraint(X);
    
    opt.box = 0;

    %constraint set
%     f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
%     X.ineq = f(opt.x);
    X = fill_constraint(X);

    opt.X = X;
    opt.scale = 0;
    out = set_path_feas(opt, order);
end

if DRAW && out.feas
    figure(1)
    clf
    syms t [1 1]
    syms x [2 1]
%     fy = f(x);
    vy = out.vval(t, x);

    hold on
    xl = [-2, 2];
    yl = [-2, 2];
    fimplicit(x'*x == Rinner^2, [xl, yl], 'k', 'DisplayName','X')
%     fimplicit(fy == 0, [xl, yl], 'DisplayName','X')

    epsilon = 1e-4;

    scatter(opt.X0(1),opt.X0(2), 100, 'ok', 'DisplayName', 'X0')
    scatter(opt.X1(1),opt.X1(2), 100, '*k', 'DisplayName', 'X1')
    
    vy0 = subs(vy, t, 0);
    fimplicit(vy0 == out.v0, [xl, yl], 'DisplayName','v(t,x) <= v(0, x0)')

    vyT = subs(vy, t, opt.Tmax);
    fimplicit(vyT == out.v1, [xl, yl], 'DisplayName','v(t, x) <= v(T, x1)')
    
%     fimplicit(x'*x == Router^2, [xl, yl], 'k', 'HandleVisibility','off')
%     fimplicit(x(1)^2 == Absx^2, [xl, yl], 'k', 'HandleVisibility','off')
    
    plot(Absx*[-1, -1], yl, 'k')

    legend('location', 'northwest')
    
    title(['Path Feasibility Certificate: time=', num2str(-out.v0, 3)])
    
    hold off
end
