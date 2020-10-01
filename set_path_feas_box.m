function [out] = set_path_feas_box(options, order)
%SET_PATH_FEAS_BOX Find a path connecting a pair of points from X0 to X1
%entirely within a set X. This is based on an optimal-control problem
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

d = 2*order;

x_yalmip = options.x;
n = length(x_yalmip);

%% Start the variables and measures

%terminal measure
mpol('tT', 1, 1);
mpol('xT', n, 1);
muT = meas([tT; xT]);

%occupation measure
mpol('t', 1, 1);
mpol('x', n, 1);
% mpol('u', n, 1);
mu = meas([t; x]);

%control and complement occupation measure
% for i = 1:n

%sigma
mpol('ts', n, 1);
mpol('xs', n, n);

%sigma_complement
mpol('tc', n, 1);
mpol('xc', n, n);

sigma = cell(n, 1);
sigma_c = cell(n, 1);
for i = 1:n
    sigma{i} = meas([ts(i); xs(:, i)]);
    sigma_c{i} = meas([tc(i); xc(:, i)]);
end
    

%assume that X0 and X1 are sets for now

%% Constraints

%support constraints
X = constraint_convert(options.X, x_yalmip, x);
if options.scale
    T = t * (1-t) >= 0;
else
    T = t * (options.Tmax - t) >= 0;
end

supp_muT = [subs(T, t, tT); xT == options.X1];
supp_mu  = [T; X];


supp_sigma = [];
supp_sigma_c = [];
for i = 1:n
    supp_sigma_curr = subs(supp_mu, [t; x], [ts(i); xs(:, i)]);
    supp_sigma_c_curr = subs(supp_mu, [t; x], [tc(i); xc(:, i)]);
    
    supp_sigma = [supp_sigma; supp_sigma_curr];
    supp_sigma_c = [supp_sigma_c; supp_sigma_c_curr];
end

%moment constraints
monT = mmon([tT; xT], d);
mon  = mmon([t; x], d);
mon_half = mmon([t; x], order);
yT = mom(monT);
% Ay = mom(diff(mon, x)*u) + mom(diff(mon, t));
       
one = ones(n, 1);

if options.scale
    scale_weight = options.Tmax;
else
    scale_weight = 1;
end

Ay_mu = mom(diff(mon, x)*-one)*scale_weight + mom(diff(mon, t));
Ay_sigma = 0;
Acont = [];

for i = 1:n
    %liouville
    mon_i = mmon([ts(i); xs(:, i)], d);
    mon_c_i = mmon([tc(i); xc(:, i)], d);
    ts_curr = ts(i);
    xs_curr = xs(:, i);

    Ay_sigma_curr = mom(2*diff(mon_i, xs(i, i)))*scale_weight;
    Ay_sigma = Ay_sigma + Ay_sigma_curr;
    
    %absolute continuity
    Acont_curr = mom(mon_i) + mom(mon_c_i) - mom(mon);
    Acont = [Acont; Acont_curr];
end

Ay = Ay_mu + Ay_sigma;

powers = genPowGlopti(n+1, d);
y0 = prod(([0; options.X0]').^powers, 2);


Liou = yT - Ay - y0;

%altogether now
%supp_con = [supp_mu0; supp_muT; supp_mu];
supp_con = [supp_muT; supp_mu; supp_sigma; supp_sigma_c];
mom_con  = [Liou == 0; Acont == 0];

objective = min(mom(tT)*scale_weight);


mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));

P = msdp(objective, ...
    mom_con, supp_con);

[status,obj_rec, m,dual_rec]= msol(P);

if options.scale
    mon_unscale = subs(mon, t, t/options.Tmax);
    mon_unscale_half = subs(mon_half, t, t/options.Tmax);
else
    mon_unscale = mon;
    mon_unscale_half = mon_half;
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
    
    
    out.Msigma = cell(n, 1);
    out.Msigma_c = cell(n, 1);
    
    y_sigma = zeros(length(out.MT), n);
    
    for i = 1:n
        out.Msigma{i} = double(mmat(sigma{i}));
        out.Msigma_c{i} = double(mmat(sigma_c{i}));
    end
    
    y_sigma = zeros(length(out.Msigma{1}), n);
    for i = 1:n
        %moments of control-occupation measure
        y_sigma(:, i) = out.Msigma{i}(1, :);                
    end
    
    %use control-occupation moments to get controllers
    u_coeff = out.Mocc \ y_sigma;
    
    u_law = u_coeff'*mon_unscale(1:length(u_coeff));
    out.u = @(ti, xi) eval(u_law, [t; x], [ti; xi]);
    out.f = @(ti, xi) eval(2*u_law - 1, [t; x], [ti; xi]);
    
    %functions    
    dual_rec_param = dual_rec{1};
    dual_rec_v = dual_rec_param(1:length(Liou));
    v = dual_rec_v'*mon_unscale;
    v1 = subs(v, x, options.X1);

    out.vval = @(ti, xi) eval(v, [t; x], [ti; xi]);    %dual v(t,x,w)
%     out.Lvval = @(ti, xi) eval(Lv, [t; x], [ti; xi]);   %Lie derivative Lv(t,x,w)
    out.v1val = @(ti) eval(v1, [t], [ti]);    %dual v(t,x,w)

    out.v0 = out.vval(0, options.X0);
    out.v1 = out.vval(options.Tmax, options.X1);
    
    %TODO: extract out the dual functions `zeta' from constraints
    
elseif status == -1
    %infeasible problem
    out.feas = 0;
end

end

