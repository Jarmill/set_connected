%Certify set disconnectedness with box control
%1d set: union of [0, low] and [high, 1]

order = 3; 
d =2*order;
T = 1; %maximum time
epsilon = 0.01; %v(0, X0) >= epsilon

%boundary of set
low = 0.5;
high = 0.7;

%initial points
% X0 = 0.2;
% X1 = 0.9;

X1 = 0.3;
X0 = 0.9;

% X0 = 0.2;
% X1 = 1;

%% variables and support sets
t = sdpvar(1,1);
x = sdpvar(1,1);

Xleft = struct('ineq', [low-x; x], 'eq', []);
Xright = struct('ineq', [1-x; x - high], 'eq', []);

X = {Xleft, Xright};
    
    
Tcons = struct('ineq', t*(T-t), 'eq', []);

All_left = struct('ineq', [Tcons.ineq; Xleft.ineq], 'eq', []);
All_right = struct('ineq', [Tcons.ineq; Xright.ineq], 'eq', []);

[v, cv] = polynomial([t; x], d);
[zeta, czeta] = polynomial([t; x], d);


%% formulate constraints
Lv = jacobian(v, t) - jacobian(v,x) - zeta;
ui = zeta + 2*jacobian(v, x);

v0 = replace(v, [t;x], [0; X0]);
vT = replace(v, [t;x], [T;X1]);



cons= [v0 >= epsilon; vT <= 0];

Tvar = 1;
vT = replace(v, [t; x], [T; X1]);

[pL, consL, coeffL] = constraint_psatz(Lv, All_left, [t;x], d);
[pR, consR, coeffR] = constraint_psatz(Lv, All_right, [t;x], d);

[pzetaL, conszetaL, coeffzetaL] = constraint_psatz(zeta, All_left, [t;x], d);
[pzetaR, conszetaR, coeffzetaR] = constraint_psatz(zeta, All_right, [t;x], d);

[puL, consuL, coeffuL] = constraint_psatz(ui, All_left, [t;x], d);
[puR, consuR, coeffuR] = constraint_psatz(ui, All_right, [t;x], d);

nonneg = [Lv; zeta; ui];


%objective

objective = 0;

%% package up
coeff = [cv; czeta; coeffL; coeffR; coeffzetaL; coeffzetaR;coeffuL; coeffuR];
cons = [cons; consL; consR; conszetaL; conszetaR; consuL; consuR];
opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);



%% plot and recovery
if sol.problem == 0

    % value(v0)

v_rec = value(cv)' * monolist([t; x], d);
zeta_rec = value(czeta)' * monolist([t; x], d);
fv = polyval_func(v_rec, [t; x]);
fzeta = polyval_func(zeta_rec, [t; x]);

    
vv0 = fv([0; X0]);
vv1star = fv([-vv0; X1]);
vv1 = fv([T; X1]);
% [vv0, vv1star, vv1]

fprintf('v(0,x0) = %0.3f, \t v(T, x1) = %0.3f \n', vv0, vv1)

figure(4)

clf
fsurf(@(t,x) fv([t;x]), [0,T,0,1], 'DisplayName', 'v(t,x)')

hold on
fcontour(@(t,x) fv([t;x]), [0,T,0,1], 'k', 'LevelList', 0, 'LineWidth', 4, 'DisplayName', 'v(t,x)=0');

scatter3(0,X0, fv([0; X0]), 400, 'ko', 'DisplayName', 'X0', 'LineWidth', 3)
scatter3(T,X1,fv([T; X1]), 400, 'k*', 'DisplayName', 'X1', 'LineWidth', 3)
% h = axes;
set(gca, 'Ydir', 'reverse')


zl = zlim;
xl_pattern = [0, T, T, 0, 0];
zl_pattern = zl([1, 1, 2, 2, 1]);
patch(xl_pattern, low*ones(5,1), zl_pattern, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'None', 'HandleVisibility', 'Off')
patch(xl_pattern, high*ones(5,1), zl_pattern, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'None', 'DisplayName', 'Region')
    
ylabel('radius')
xlabel('time')
zlabel('v')
title(sprintf('Auxiliary Function on 1d Separation (order=%d, box)', order), 'FontSize', 16)
legend('location', 'northwest')
end