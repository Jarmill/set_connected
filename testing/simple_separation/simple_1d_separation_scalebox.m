%Certify set disconnectedness
% %union of [0, low] and [high, 1]

order =5; 
d =2*order;
T = 1; %maximum time
% r = 5;
r = 0;
% alpha = 10;
alpha = 0;
% epsilon = 0.01;
epsilon = 1;

xbox = [-1; 1];
% xbox = [-1; 0.5];
% xbox = [-0.3; 0.3];
% xbox = [-2;2];

% low = 0.4;
% low = 0.7;
% high = 0.8;
% low_alpha= 0.5;
% high_alpha= 0.55;
low_alpha = 0.7;
high_alpha = 0.75;
low = xbox(2)*low_alpha+ xbox(1)*(1-low_alpha);
high = xbox(2)*high_alpha+ xbox(1)*(1-high_alpha);



% low = 0.6;
% high = 0.8;
% 
% X0 = 0.3;
% X1 = 0.9;

% X0_alpha = 0.2;
% X0_alpha = 0.3;
X0_alpha = 0.35;
X1_alpha = 0.9;
% X1_alpha = 0.85;

X0 = xbox(2)*X0_alpha + xbox(1)*(1-X0_alpha);
X1 = xbox(2)*X1_alpha + xbox(1)*(1-X1_alpha);

%% variables and support sets
t = sdpvar(1,1);
x = sdpvar(1,1);
u = sdpvar(1,1);

Xleft = struct('ineq', [(low-x)*(x-xbox(1))], 'eq', []);
Xright = struct('ineq', [(xbox(2)-x)*(x - high)], 'eq', []);

X = {Xleft, Xright};
    
Tcons = struct('ineq', t*(1-t), 'eq', []);

% U = struct('ineq', [1-u; u+1], 'eq', []);
U = struct('ineq', 1-u^2, 'eq', []);

All_left = struct('ineq', [Tcons.ineq; Xleft.ineq; U.ineq], 'eq', []);
All_right = struct('ineq', [Tcons.ineq; Xright.ineq; U.ineq], 'eq', []);

[v, cv] = polynomial([t; x], d);

%test sos program 

%% formulate constraints
Lv = jacobian(v, t) + jacobian(v,x)*u*T + alpha*v;

v0 = replace(v, [t;x], [0; X0]);
vT = replace(v, [t;x], [1;X1]);

cons= [v0 >= epsilon; vT <= 0];

% % vT = replace(v, [t; x], [T; X1]);

% Rx = (1+t+x)^r;
Ru = (1+t+x+u^2)^r;

[pL, consL, coeffL] = constraint_psatz(Ru*Lv, All_left, [t; x; u], d);
[pR, consR, coeffR] = constraint_psatz(Ru*Lv, All_right, [t; x; u], d);

%objective

objective = norm(cv);

%% package up
coeff = [cv; coeffL; coeffR];
cons = [cons; consL; consR];
opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);



% value(v0)
%% plot and recovery


if sol.problem == 0
    v_rec = value(cv)' * monolist([t; x], d);
    nv_rec = norm(value(cv));
fv = polyval_func(v_rec, [t; x]);

vv0 = fv([0; X0]);
vv1star = fv([-vv0; X1]);
vv1 = fv([1; X1]);
% [vv0, vv1star, vv1]
fprintf('v(0,x0) = %0.3f, \t v(T, x1) = %0.3f, \t norm(cv) = %0.3f \n', vv0, vv1, nv_rec)

figure(5)

clf
fsurf(@(t,x) fv([t;x]), [0,1,xbox'], 'DisplayName', 'v(t,x)')

hold on
fcontour(@(t,x) fv([t;x]), [0,1,xbox'], 'k', 'LevelList', 0, 'LineWidth', 4, 'DisplayName', 'v(t,x)=0');

scatter3(0,X0, fv([0; X0]), 400, 'ko', 'DisplayName', 'X0', 'LineWidth', 3)
scatter3(1,X1,fv([1; X1]), 400, 'k*', 'DisplayName', 'X1', 'LineWidth', 3)
% h = axes;
set(gca, 'Ydir', 'reverse')


zl = zlim;
xl_pattern = [0, 1, 1, 0, 0];
zl_pattern = zl([1, 1, 2, 2, 1]);
patch(xl_pattern, low*ones(5,1), zl_pattern, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'None', 'HandleVisibility', 'Off')
patch(xl_pattern, high*ones(5,1), zl_pattern, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'None', 'DisplayName', 'Region')
    
ylabel('radius')
xlabel('time')
zlabel('v')
title(sprintf('Barrier Function on 1d Separation (order=%d)', order), 'FontSize', 16)
legend('location', 'northwest')
    disp('feasible')
else
    disp('infeasible')
end