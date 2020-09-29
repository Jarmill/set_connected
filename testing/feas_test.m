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

    FEAS = 0;
    SET = 0;
    if FEAS
%         opt.X0 = [1.25; -1];
%         opt.X1 = [1.5; 0.5];
        
         opt.X0 = [1.05; -1];
         opt.X1 = [0.4; 0.7];
        
        
%         opt.X0 = [1.5; -1];
%         opt.X1 = [1.6; 1]; %x=1.5 works, 1.6 does not
    else
            opt.X0 = [-0.75; 0];
            opt.X1 = [1.5; 0.5];
%             opt.X1 = [1; 0.5];
%             opt.X1 = [1; -0.5];
    end

    opt.box = 0;

    %constraint set
    f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
    X.ineq = f(opt.x);
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
    fy = f(x);
    vy = out.vval(t, x);

    hold on
    xl = [-2, 2];
    yl = [-2, 2];
    fimplicit(fy == 0, [xl, yl], 'DisplayName','X')

    epsilon = 1e-4;

    scatter(opt.X0(1),opt.X0(2), 100, 'ok', 'DisplayName', 'X0')
    scatter(opt.X1(1),opt.X1(2), 100, '*k', 'DisplayName', 'X1')
    
    vy0 = subs(vy, t, 0);
    fimplicit(vy0 == out.v0, [xl, yl], 'DisplayName','v(t,x) <= v(0, x0)')

    vyT = subs(vy, t, opt.Tmax);
    fimplicit(vyT == out.v1, [xl, yl], 'DisplayName','v(t, x) <= v(T, x1)')
    
    legend('location', 'northwest')
    
    title(['Path Feasibility Certificate: time=', num2str(-out.v0, 3)])
    
    hold off
end
