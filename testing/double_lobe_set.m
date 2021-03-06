%try to do proofs of disconnectedness of double-lobe region

n = 2;

t = sdpvar(1,1);
x = sdpvar(n,1);
u = sdpvar(n,1);


order = 1;
d =2*order;
T = 2; %maximum time

FEAS = 0;

% if FEAS
%     x0 = [1.25; -1];
%     x1 = [1.5; 0.5];
% else
    x0 = [-0.75; -0.25];
    x1 = [1.5; 0.5];
% end

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

Tvar = 1;
vT = replace(v, [t; x], [T; x1]);
if Tvar
    [pT, consT, coeffT] = constraint_psatz(v1t + t, Tcons, t, d);
else
    pT = vT;
    coeffT = [];
    consT = (pT >= 0);
end
[pL, consL, coeffL] = constraint_psatz(-Lv, Allcons, [t; x; u], d);

%objective

% objective= -v0;
objective = v0;

%package up
coeff = [cv; coeffT; coeffL];
cons = [consT; consL];
opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;
if FEAS
    [sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);
else
    cons = [cons; v0 == 1];
    [sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);
end
value(v0)