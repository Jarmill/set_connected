function [out] = set_path_infeas(options, order)
%SET_PATH_INFEAS Check for infeasibility certificate of path connection
%Use Farkas Lemma
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;


%try to do proofs of disconnectedness of double-lobe region

%This is a farkas certificate of infeasibility
%there exists a function v where v <= 0 on X0, v >= 0 on X1, and v
%decreases along trajectories. Trajectories will never cross the surface
%v(t,x) = 0, and this level sets separates X0 and X1.


%% Process options and variables

%variables of problems
x = options.x;
n = length(x);
u = sdpvar(n, 1);
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

U = struct;
if strcmp(options.U, 'box')
    U.ineq = 1 - u.^2;
else
    U.ineq = 1 - u'*u;
end
U = fill_constraint(U);

Allcons.ineq = [T.ineq; X.ineq; U.ineq];
Allcons.eq = X.eq;
Allcons = fill_constraint(Allcons);


%% form constraints in Yalmip
[v, cv] = polynomial([t; x], d);

%test sos program 

%formulate constraints
% Lv = jacobian(v, [t; x])*[1; u];


%value function at (0, X0) is negative
if isstruct(X0)    
    %X0 is a set
    set0 = true;
    X0 = fill_constraint(X0);
    v0 = replace(v, t, 0);
    [p0, cons0, coeff0] = constraint_psatz(-v0 - 1, X0, x, d);
else
    %X0 is a point    
    set0 = false;
    v0 = replace(v, [t;x], [0; X0]);
    cons0 = (v0 == -1);
    coeff0 = [];    
end



%value function on (T, X1) is positive
if isstruct(X1)
    %X1 is a set
    set1 = true;
    X1 = fill_constraint(X1);
    if options.scale
        vT = replace(v, t, t);
    else
        vT = replace(v, t, options.Tmax);
    end
    
    [pT, consT, coeffT] = constraint_psatz(vT, X1, x, d);
else
    %X1 is a point    
    set1 = false;
    if options.scale
        vT = replace(v, [t; x], [1; X1]);
    else
        vT = replace(v, [t; x], [options.Tmax; X1]);
    end
%     [pT, consT, coeffT] = constraint_psatz(vT, X1, x, d);
    consT = (vT >= 0);
    coeffT = [];
end

%value function v decreases along trajectories
%therefore impossible for trajectories to go from X0 to X1
if options.scale
    Lv = jacobian(v, t) + options.Tmax * jacobian(v,x)*u;
else    
    Lv = jacobian(v, t) + jacobian(v,x)*u;
end
[pL, consL, coeffL] = constraint_psatz(-Lv, Allcons, [t; x; u], d);




objective = 0;

%% form the optimization problem
coeff = [cv; coeff0; coeffT; coeffL];
cons = [cons0; consT; consL];
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
    Lv_rec = jacobian(v_rec, t) + jacobian(v_rec,x)*u;
    vval = polyval_func(v_rec, [t; x]);
    Lvval = polyval_func(Lv_rec, [t; x; u]);
    
    out.vval = vval;
    out.Lvval = Lvval;
    
    if set0       
        out.v0 = value(v0);
    end
    if set1
        out.vT = value(vT);
    end
    
elseif sol.problem == 1
    %infeasible, there may be a path from X0 to X1
    out.farkas = false;
end



end

