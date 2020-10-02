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


    order = 3;
    d =2*order;
    T = 3; %maximum time

    FEAS = 0;

    
    if FEAS
        opt.X0 = [1; -1];
        opt.X1 = [1; 1];
    else
        opt.X0 = [1; -1];
        opt.X1 = [-1; -1];
    end

    %constraint set
%     Rinner = 1.3;
%     Rinner = 1;
    Router = 1.8;
    Absx = 0.5;

    x = opt.x;
    X = struct;
%     if ~SEP
        %annulus is connected
%         X.ineq = [x'*x - Rinner^2; Router^2 - x'*x];
%     else
        %disk is split
        X.ineq = [Router^2 - x'*x; x(1)^2 - Absx^2];
%     end
    X = fill_constraint(X);
    
    opt.box = 0;

    %constraint set
%     f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
%     X.ineq = f(opt.x);
    X = fill_constraint(X);

    opt.X = X;
    opt.scale = 0;
    
    BOX = 1;
    
    if BOX
        out = set_path_infeas_box(opt, order);
    else
        out = set_path_infeas(opt, order);
    end
      
end

if DRAW && out.farkas
    figure(1)
    clf
    syms t [1 1]
    syms x [2 1]
%     fy = f(x);
    vy = out.vval([t; x]);

    hold on
    xl = [-2, 2];
    yl = [-2, 2];
    
    % TODO: faster plots
    %explicitly plot the circles and lines of the annulus
    %make it easy
    N = 200;
    theta = linspace(0, 2*pi, N);
    circ = [cos(theta); sin(theta)];
    
%     plot(Rinner*circ(1, :), Rinner*circ(2, :), 'k', 'DisplayName','X')
%     fimplicit(x'*x == Rinner^2, [xl, yl], 'k', 'DisplayName','X')
%     fimplicit(fy == 0, [xl, yl], 'DisplayName','X')

    epsilon = 1e-4;

    scatter(opt.X0(1),opt.X0(2), 100, 'ok', 'DisplayName', 'X0')
    scatter(opt.X1(1),opt.X1(2), 100, '*k', 'DisplayName', 'X1')
    
%     vy0 = subs(vy, t, 0);
%     fimplicit(vy0 == out.v0, [xl, yl], 'DisplayName','v(t,x) <= v(0, x0)')
% 
%     vyT = subs(vy, t, opt.Tmax);
%     fimplicit(vyT == out.v1, [xl, yl], 'DisplayName','v(t, x) <= v(T, x1)')
%     

    epsilon = 1e-4;
    vy0 = subs(vy, t, 0);
    color0 = [0.4940, 0.1840, 0.5560];
    colorT = [0.4660, 0.6740, 0.1880];


    fimplicit(vy0 == -1, [xl, yl], 'DisplayName','v(0, x) <= -1' ,'Color', color0)
    fimplicit(vy0 == -epsilon, [xl, yl], ':', 'Color', color0, 'DisplayName','v(0, x) < 0')    


    vyT = subs(vy, t, opt.Tmax);
    fimplicit(vyT == 1, [xl, yl], 'DisplayName','v(T, x) >= 1', 'Color', colorT)
    fimplicit(vyT == epsilon, [xl, yl], ':', 'Color', colorT, 'DisplayName','v(T, x) > 0')

    
%     fimplicit(x'*x == Router^2, [xl, yl], 'k', 'HandleVisibility','off')\
    



    plot(Router*circ(1, :), Router*circ(2, :), 'k', 'DisplayName','X')
    plot(-Absx*[1,1], yl, 'k', 'HandleVisibility','off')
    plot(Absx*[1,1], yl, 'k', 'HandleVisibility','off')
    
%     fimplicit(x(1)^2 == Absx^2, [xl, yl], 'k', 'HandleVisibility','off')
    
    

    legend('location', 'northwest')
    
    title(['Farkas Infeasibility Certificate')
    
    hold off
end
