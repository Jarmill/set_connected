x = sdpvar(2,1);

order = 4;
d =2*order;


%Points in left lobe
%[-1; 0]
%[-0.75; -0.25]
%[-1.25; 0.25]

%Points on right lobe
%[1.25; -1]
%[1; 1]
%[1.5; 0.5]

%constraint set
f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
X.ineq = f(x);
X = fill_constraint(X);

%obj = x(1);
obj = x(2);

%test sos program 
%minimize obj on g(x) >= 0
gamma = sdpvar(1,1);

% [v, cv] = polynomial(vars, d);
[p0, cons, coeff] = constraint_psatz(obj - gamma, X, x, d);


opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;

objective = -gamma;
[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff; gamma]);

peak_val = -value(gamma);



syms xv [2 1];
fv = f(xv);
% fimplicit(fv);