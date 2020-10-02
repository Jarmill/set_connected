function [out] = set_path_infeas_box_indep(options, order)
%SET_PATH_INFEAS_BOX_INDEP Check for infeasibility certificate of path connection
%Use Farkas Lemma to find a violating polynomial
%
%Efficient formulation where U = [-1, 1]^n
%unsigned model, U warped to [0, 1]^n
%time independent formulation
%
%Input:
%   options: set_path_options data structure
%   order:   order of relaxation (degree = 2*order)
%
%Output:
%   out:     out data structure
%       farkas: whether a farkas certificate has been found (1 or 0)
%       vval:   A function that evaluates v(t,x)
%       Lvval:  A function that evaluates Luv = dv/dt + u'*dv/dx
%       v0:     If X0 is a set, v0 is the value of v(0, x0)
%       v1:     If X1 is a set, v0 is the value of v(T, x1)

%This is a farkas certificate of infeasibility
%there exists a function v where v <= 0 on X0, v >= 0 on X1, and v
%decreases along trajectories. Trajectories will never cross the surface
%v(t,x) = 0, and this level sets separates X0 and X1.


%% Process options and variables
%variables of problems
x = options.x;
n = length(x);



d =2*order;
% T = options.Tmax; %maximum time

%constraints
X0 = options.X0;
X1 = options.X1;

X = options.X; 

%time

Allcons.ineq = [X.ineq];
Allcons.eq = X.eq;
Allcons = fill_constraint(Allcons);


%% form constraints in Yalmip
[v, cv] = polynomial(x, d);

zeta = [];
coeff_zeta = [];
cons_zeta_sos = [];
zeta_sum = 0;
for i = 1:n
    [pzeta, czeta] = polynomial(x, d);
    coeff_zeta = [coeff_zeta; czeta];
    cons_zeta_sos = [cons_zeta_sos; sos(pzeta)];
    
    zeta = [zeta; pzeta];
    zeta_sum = zeta_sum + pzeta;
end

%test sos program 

%formulate constraints
% Lv = jacobian(v, x)*[1; u];


%value function at (0, X0) is negative
v0 = replace(v, x, X0);
vT = replace(v, x, X1);

cons0 = (v0 == -1);
consT = (vT  == 1);

%value function v decreases along trajectories
%therefore impossible for trajectories to go from X0 to X1


%revise this for the box code


div_v = jacobian(v,x)*ones(n, 1);
Lv = -div_v - zeta_sum;

[pL, consL, coeffL] = constraint_psatz(-Lv, Allcons, x, d);

cons_u = [];
coeff_u = [];
for i = 1:n
    term_ui = zeta(i) + 2*jacobian(v,x(i));
    [p_ui, cons_ui, coeff_ui] = constraint_psatz(term_ui, Allcons, x, d);
    
    cons_u = [cons_u; cons_ui];
    coeff_u = [coeff_u; coeff_ui];
end

objective = 0;
% objective = -cv(1);

%% form the optimization problem
coeff = [cv; coeffL; coeff_u; coeff_zeta];
cons = [cons0; consT; consL; cons_zeta_sos; cons_u];
opts = sdpsettings('solver', options.solver);
opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);

%% Evaluate and package output
out = struct;


mlist = monolist(x, d);


%test for strict feasibility
min_eig = cellfun(@(c) min(eig(c)), Gram);
if min(min_eig) < 0
    %infeasible
    sol.problem = 1;
end

if sol.problem == 0
    %feasible, there exists a farkas certificate such that X0 and X1 are
    %separated in X
    out.farkas = true;
    
    %find the value function v(t,x) 
    v_rec = value(cv)' * mlist;
%     Lv_rec = jacobian(v_rec, t) + jacobian(v_rec,x)*u;
    vval = polyval_func(v_rec, x);
%     Lvval = polyval_func(Lv_rec, [t; x; u]);
    zetaval = 0;
    
    coeff_zeta_val = reshape(value(coeff_zeta),[], 2);
    zeta_rec = coeff_zeta_val'*mlist;
    zetaval = polyval_func(zeta_rec, x);

    out.vval = @(xi) vval(xi(2:end));
    out.zetaval = @(xi) zetaval(xi(2:end));

    out.v0 = value(v0);
    out.vT = value(vT);
    
    out.min_eig = min_eig;
elseif sol.problem == 1
    %infeasible, there may be a path from X0 to X1
    out.farkas = false;
    out.min_eig = min_eig;
end



end

