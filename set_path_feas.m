function [out] = set_path_feas(options, order)
%SET_PATH_FEAS Find a path connecting a pair of points from X0 to X1
%entirely within a set X. This is based on an optimal-control problem
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

d = 2*order;

x_yalmip = options.x;
n = length(x_yalmip);

%% Start the variables and measures
mpol('tT', 1, 1);
mpol('xT', n, 1);
muT = meas([tT; xT]);

mpol('t', 1, 1);
mpol('x', n, 1);
mpol('u', n, 1);
mu = meas([t; x; u]);


%assume that X0 and X1 are sets for now

%% Constraints

%support constraints
X = constraint_convert(options.X, x_yalmip, x);

if strcmp(options.U, 'box')
    U = 1 - u.^2 >= 0;
else
    U= 1 - u'*u >= 0;
end

if options.scale
    supp_muT = [tT*(1-tT) >= 0; xT == options.X1];
    supp_mu = [t*(1-t) >= 0; X; U];
else
    supp_muT = [tT*(options.Tmax-tT) >= 0; xT == options.X1];
    supp_mu = [t*(options.Tmax-t) >= 0; X; U];
end

%moment constraints
monT = mmon([tT; xT], d);
mon  = mmon([t; x], d);

yT = mom(monT);
Ay = mom(diff(mon, x)*u) + mom(diff(mon, t));
       

powers = genPowGlopti(n+1, d);
y0 = prod(([0; options.X0]').^powers, 2);


Liou = yT - Ay - y0;

%altogether now
%supp_con = [supp_mu0; supp_muT; supp_mu];
supp_con = [supp_muT; supp_mu];
mom_con  = (Liou == 0);

objective = min(mom(tT));


mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));

P = msdp(objective, ...
    mom_con, supp_con);

[status,obj_rec, m,dual_rec]= msol(P);

if options.scale
    mon_unscale = subs(mon, t, t/options.Tmax);
else
    mon_unscale = mon;
end

out = struct;
if status == 0
    %feasible problem
    out.feas = 1;
    
    %moment matrices
    powers_half = genPowGlopti(n+1, order);
    y0_half = prod(([0; options.X0]').^powers_half, 2);
    out.M0 = y0_half*y0_half';


    out.MT = double(mmat(muT));
    out.Mocc = double(mmat(mu));
    
    %functions    
    dual_rec_v = dual_rec{1};
    v = dual_rec_v'*mon_unscale;
    v1 = subs(v, x, options.X1);

    out.vval = @(ti, xi) eval(v, [t; x], [ti; xi]);    %dual v(t,x,w)
    out.Lvval = @(ti, xi) eval(Lv, [t; x], [ti; xi]);   %Lie derivative Lv(t,x,w)
    out.v1val = @(ti) eval(v1, [t], [ti]);    %dual v(t,x,w)

    out.v0 = out.vval(0, options.X0);
    out.v1 = out.vval(options.Tmax, options.X1);
elseif status == -1
    %infeasible problem
    out.feas = 0;
end

end

