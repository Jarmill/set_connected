function [out] = set_path_infeas_box(options, order)
%SET_PATH_INFEAS_BOX Check for infeasibility certificate of path connection
%Use Farkas Lemma to find a violating polynomial
%
%Efficient formulation where U = [-1, 1]^n
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
t = sdpvar(1, 1);




d =2*order;
% T = options.Tmax; %maximum time

%constraints
X0 = options.X0;
X1 = options.X1;

X = options.X; 

T = struct;
if options.scale
    T.ineq = t*(1-t);
else
    T.ineq = t*(options.Tmax-t);
end
T = fill_constraint(T);

Allcons.ineq = [T.ineq; X.ineq];
Allcons.eq = X.eq;
Allcons = fill_constraint(Allcons);


%% form constraints in Yalmip
[v, cv] = polynomial([t; x], d);

zeta = [];
coeff_zeta = [];
cons_zeta_sos = [];
zeta_sum = 0;
for i = 1:n
    [pzeta, czeta] = polynomial([t; x], d);
    coeff_zeta = [coeff_zeta; czeta];
    cons_zeta_sos = [cons_zeta_sos; sos(pzeta)];
    
    zeta = [zeta; pzeta];
    zeta_sum = zeta_sum + pzeta;
end

%test sos program 

%formulate constraints
% Lv = jacobian(v, [t; x])*[1; u];


%value function at (0, X0) is negative
% if isstruct(X0)    
%     %X0 is a set
%     set0 = true;
%     X0 = fill_constraint(X0);
%     v0 = replace(v, t, 0);
%     [p0, cons0, coeff0] = constraint_psatz(-v0 - 1, X0, x, d);
% else
    %X0 is a point    
    set0 = false;
    v0 = replace(v, [t;x], [0; X0]);
%     cons0 = (v0 == 1);
    cons0 = (v0 -1 >= 0);
    coeff0 = [];    
% end



%value function on (T, X1) is positive
% if isstruct(X1)
%     %X1 is a set
%     set1 = true;
%     X1 = fill_constraint(X1);
%     if options.scale
%         vT = replace(v, t, t);
%     else
%         vT = replace(v, t, options.Tmax);
%     end
%     
%     [pT, consT, coeffT] = constraint_psatz(vT, X1, x, d);        
% else
    %X1 is a point    
    set1 = false;
    if options.scale
        vT = replace(v, [t; x], [1; X1]);
    else
        vT = replace(v, [t; x], [options.Tmax; X1]);
    end
%     [pT, consT, coeffT] = constraint_psatz(vT, X1, x, d);
%     consT = (vT == -1);
    consT = -vT -1>= 0; 
    coeffT = [];
% end

%value function v decreases along trajectories
%therefore impossible for trajectories to go from X0 to X1


%revise this for the box code

if options.scale
    scale_weight = options.Tmax;
else
    scale_weight = 1;
end


Lv = jacobian(v, t) - scale_weight*jacobian(v,x)*ones(n, 1) - zeta_sum;

% [pL, consL, coeffL] = constraint_psatz(-Lv, Allcons, [t; x], d);

[pL, consL, coeffL] = constraint_psatz(Lv, Allcons, [t; x], d);

cons_u = [];
coeff_u = [];
for i = 1:n
    term_ui = zeta(i) + scale_weight*2*jacobian(v,x(i));
    [p_ui, cons_ui, coeff_ui] = constraint_psatz(term_ui, Allcons, [t; x], d);
    
    cons_u = [cons_u; cons_ui];
    coeff_u = [coeff_u; coeff_ui];
end

objective = 0;

%% form the optimization problem
coeff = [cv; coeff0; coeffT; coeffL; coeff_u; coeff_zeta];
cons = [cons0; consT; consL; cons_zeta_sos; cons_u];
opts = sdpsettings('solver', options.solver);
opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);

%% Evaluate and package output
out = struct;

if options.scale
    mlist = monolist([t/options.Tmax; x], d);
else
    mlist = monolist([t; x], d);
end


if sol.problem == 0
    %feasible, there exists a farkas certificate such that X0 and X1 are
    %separated in X
    out.farkas = true;
    
    %find the value function v(t,x) 
    v_rec = value(cv)' * mlist;
%     Lv_rec = jacobian(v_rec, t) + jacobian(v_rec,x)*u;
    vval = polyval_func(v_rec, [t; x]);
%     Lvval = polyval_func(Lv_rec, [t; x; u]);
    zetaval = 0;
    
    coeff_zeta_val = reshape(value(coeff_zeta),[], 2);
    zeta_rec = coeff_zeta_val'*mlist;
    zetaval = polyval_func(zeta_rec, [t; x]);
    
    out.vval = vval;
    out.zetaval = zetaval;
%     out.Lvval = Lvval;
    
    if ~set0       
        out.v0 = value(v0);
    end
    if ~set1
        out.v1 = value(vT);
    end
    
else
    %sol.problem == 1   infeasible
    %sol.problem == 4   numerical problems
    
    %infeasible, there may be a path from X0 to X1
    out.farkas = false;
end



end

