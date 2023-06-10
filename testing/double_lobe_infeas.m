%try to do proofs of disconnectedness of double-lobe region

%This is a farkas certificate of infeasibility
%there exists a function v where v <= 0 on X0, v >= 0 on X1, and v
%decreases along trajectories. Trajectories will never cross the surface
%v(t,x) = 0, and this level sets separates X0 and X1.


n = 2;

t = sdpvar(1,1);
x = sdpvar(n,1);
u = sdpvar(n,1);


order = 1;
d =2*order;
T = 2; %maximum time

FEAS = 0;

if FEAS
    x0 = [1.25; -1];
    x1 = [1.5; 0.5];
else
    x0 = [-0.75; -0.25];
    x1 = [1.5; 0.5];
end

%constraint set
f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
X.ineq = f(x);
X = fill_constraint(X);

Tcons.ineq = t*(T-t);
Tcons = fill_constraint(Tcons);

U.ineq = 1 - u'*u;

Allcons.ineq = [Tcons.ineq; X.ineq; U.ineq];
Allcons = fill_constraint(Allcons);

[v, cv] = polynomial([t; x], d);

%test sos program 

%formulate constraints
% Lv = jacobian(v, [t; x])*[1; u];
Lv = jacobian(v, t) + jacobian(v,x)*u;

v0 = replace(v, [t;x], [0; x0]);
v1t = replace(v, x, x1);

cons= [];

Tvar = 0;
vT = replace(v, [t; x], [T; x1]);

pT = vT;
coeffT = [];
consT = (pT >= 0);

[pL, consL, coeffL] = constraint_psatz(-Lv, Allcons, [t; x; u], d);

%objective

% objective= -v0;
objective = v0;

%package up
coeff = [cv; coeffT; coeffL];
cons = [consT; consL];
opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;

cons = [cons; v0 == -1];
[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);

% value(v0)




v_rec = value(cv)' * monolist([t; x], d);
fv = polyval_func(v_rec, [t; x]);

vv0 = fv([0; x0]);
vv1 = fv([T; x1]);
[vv0, vv1];