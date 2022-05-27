%Certify set disconnectedness
% %union of [0, low] and [high, 1]

order = 1; 
d =2*order;
T = 1; %maximum time
% r = 5;
r = 0;
% epsilon = 0.01;
epsilon = 1;


xbox = sdpvar(1,2);
% xbox = [-1; 1];
% xbox = [-1; 0.5];
% xbox = [-0.3; 0.3];
% xbox = [-2;2];

% low = 0.4;
% low = 0.7;
% high = 0.8;
% low_alpha= 0.5;
% high_alpha= 0.55;
low_alpha = 0.7;
high_alpha = 0.71;
low = xbox(2)*low_alpha+ xbox(1)*(1-low_alpha);
high = xbox(2)*high_alpha+ xbox(1)*(1-high_alpha);



% low = 0.6;
% high = 0.8;
% 
% X0 = 0.3;
% X1 = 0.9;

X0_alpha = 0.2;
X1_alpha = 0.9;

X0 = xbox(2)*X0_alpha + xbox(1)*(1-X0_alpha);
X1 = xbox(2)*X1_alpha + xbox(1)*(1-X1_alpha);

%% variables and support sets
t = sdpvar(1,1);
x = sdpvar(1,1);
u = sdpvar(1,1);

Xleft = struct('ineq', [(low-x)*(x-xbox(1))], 'eq', []);
Xright = struct('ineq', [(xbox(2)-x)*(x - high)], 'eq', []);

X = {Xleft, Xright};
    
Tcons = struct('ineq', t*(T-t), 'eq', []);

% U = struct('ineq', [1-u; u+1], 'eq', []);
U = struct('ineq', 1-u^2, 'eq', []);

All_left = struct('ineq', [Tcons.ineq; Xleft.ineq; U.ineq], 'eq', []);
All_right = struct('ineq', [Tcons.ineq; Xright.ineq; U.ineq], 'eq', []);

[v, cv] = polynomial([t; x], d);

%test sos program 

%% formulate constraints
Lv = jacobian(v, t) + jacobian(v,x)*u;

v0 = replace(v, [t;x], [0; X0]);
vT = replace(v, [t;x], [T;X1]);

cons= [(v0 >= epsilon):'initial'; (vT <= 0):'terminal'];
% % vT = replace(v, [t; x], [T; X1]);


[pL, consL, GramL, MuL] = psatz(Lv, All_left, order, [t;x;u]);
[pR, consR, GramR, MuR] = psatz(Lv, All_right, order, [t;x;u]);
% [pL, consL, GramL, MuL] = psatz(Lv, All_left, order, [t;x;u])

% [pL, consL, coeffL] = constraint_psatz(Ru*Lv, All_left, [t; x; u], d);
% [pR, consR, coeffR] = constraint_psatz(Ru*Lv, All_right, [t; x; u], d);

%objective

objective = norm(cv);

%% package up
% coeff = [cv; coeffL; coeffR];
cons = [cons; consL; consR];
opts = sdpsettings('solver', 'mosek');

P = optimizer(cons,objective,opts,xbox, cv);

cv_rec = P([-1,1]);
% P([-0.3,0.3])

% sol = optimize(cons, objective, opts);
% opts.sos.model = 2;

% value(v0)
%% plot and recovery


if sol.problem == 0
    v_rec = value(cv)' * monolist([t; x], d);
    nv_rec = norm(value(cv));
fv = polyval_func(v_rec, [t; x]);

vv0 = fv([0; X0]);
vv1star = fv([-vv0; X1]);
vv1 = fv([T; X1]);
% [vv0, vv1star, vv1]
fprintf('v(0,x0) = %0.3f, \t v(T, x1) = %0.3f, \t norm(cv) = %0.3f \n', vv0, vv1, nv_rec)

% nGram = size(GramL,1);
% Gram_rec = cell(nGram, 2);
% for i = 1:nGram
%     Gram_rec{i,1} = value(GramL{i, 1});
%     Gram_rec{i,2} = value(GramR{i, 1});
% end
figure(5)

clf
fsurf(@(t,x) fv([t;x]), [0,T,xbox'], 'DisplayName', 'v(t,x)')

hold on
fcontour(@(t,x) fv([t;x]), [0,T,xbox'], 'k', 'LevelList', 0, 'LineWidth', 4, 'DisplayName', 'v(t,x)=0');

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
title(sprintf('Auxiliary Function on 1d Separation (order=%d)', order), 'FontSize', 16)
legend('location', 'northwest')
else
    fprintf('infeasible\n')
end