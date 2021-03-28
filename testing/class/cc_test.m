SOLVE = 1;
PLOT = 1;


opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
opt.Tmax = 2;


opt.scale = 0;

order = 3;
d = 2*order;

% X0_feas = [1.25; -1];
% X1_feas = [1.5; 0.5];

X_bump =full(sparse([1 2], [1 3], [-0.2, -0.2]));

% X0_infeas = [-0.75; 0.5];
X0_infeas = [-0.75; 0.5]+ X_bump;

%points
X1_infeas = [1.5; 0.5] + X_bump;

%small circle
% X1_infeas.ineq  = 0.2^2 - sum((opt.x - [1.5; 0.5]).^2);

f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
X.ineq = f(opt.x);
X = fill_constraint(X);

opt.X = X;
opt.X0 = X0_infeas;
opt.X1 = X1_infeas;

if SOLVE
IM = infeas_manager(opt);
out = IM.infeas(d);
end

if PLOT && out.farkas
    
    figure(1)
    clf
    syms tv [1 1];
    syms xv [2 1];
    fx = f(xv);
    v0x = out.func.v0(xv);
    v1x = out.func.v1(xv);
    hold on
    xl = [-3, 2.5];
    yl = [-2, 2];
    if opt.scale
        tl = [0, 1];
    else
        tl = [0, opt.Tmax];
    end
    
    color0 = [0.4940, 0.1840, 0.5560];
    colorT = [0.4660, 0.6740, 0.1880];
    
%     for i = 1:size(X0_infeas, 1)
    scatter(opt.X0(1, :), opt.X0(2, :), 100, 'ok', 'DisplayName', 'X0')
    scatter(opt.X1(1, :), opt.X1(2, :), 100, '*k', 'DisplayName', 'X1')
    
    fimplicit(v0x == 0, [xl, yl],':', 'DisplayName','v(0, x) = 0', 'Color', color0)
    fimplicit(v1x == 0, [xl, yl],':', 'DisplayName','v(T, x) = 0', 'Color', colorT)
    
    fimplicit(v0x == 1, [xl, yl],'DisplayName','v(0, x) = 1', 'Color', color0)
    fimplicit(v1x == -1, [xl, yl],'DisplayName','v(T, x) = -1', 'Color', colorT)
    
    fimplicit(fx, [xl, yl], 'k', 'DisplayName','X')
        
    legend('location', 'northwest')
end
% [cc_all, poly_var] = IM.make_program(d);


% %TODO: 
% %write solving code in class
% opts = sdpsettings('solver', opt.solver);
% opts.sos.model = 2;
% [sol, monom, Gram, residual] = solvesos(cc_all.con, 0, opts, cc_all.coef);
% 
% 
% %recovery code in class
% %evaluation of functions
% %visualization
% [cv,mv] = coefficients(poly_var.v,[poly_var.t; poly_var.x]);
% v_eval = value(cv)'*mv;
% 
% [cz, mz] = coefficients(poly_var.zeta,[poly_var.t; poly_var.x]);
% zeta_eval = value(cz)*mz;