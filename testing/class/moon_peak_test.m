%moon test 
%there was a bug in the code, X was set to {} before adding circles
%attempt to use peak estimation instead of set connectedness


SOLVE = 1;
PLOT = 1;
SAMPLE = 0;
EVAL = 0;
FEAS = 0;

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
t = opt.t;
x = opt.x;
opt.Tmax = 2;


opt.scale = 1;
opt.verbose = 1;

%moon set
inner_rad = 0.7;
inner_x = 0.35;

X_moon = struct;
X_moon.ineq = [1 - opt.x(1)^2 - opt.x(2)^2; 
                (opt.x(1)-inner_x)^2 + opt.x(2)^2 - inner_rad^2];
% X_moon = fill_constraint(X_moon);

X = {X_moon};

%circles
MULTIPLE_CIRCLES = 0;

if MULTIPLE_CIRCLES
    circle_rad = [0.4; 0.3; 0.3];
    circle_center = [0.4, 0;
                 -1, 1;
                 -1, -1]';

else    
    circle_rad = [0.4];
    circle_center = [0.4, 0]';
end


X_circ = {};
for i = 1:length(circle_rad)
%     X_circ = struct;
    rx = circle_rad(i);
    cx = circle_center(:, i);
    X_circ.ineq = rx^2 - (opt.x(1)-cx(1))^2 - (opt.x(2)-cx(2))^2;
    X = [X; {X_circ}];
end


% X0_feas = ;
% X1_feas = [0; -0.9];

% X0_infeas = [0 0 0.7 0.7 -0.8; 0.8 -0.8 0.65 -0.65 0 ];
X0_infeas = [0; 0.9];
X1_infeas = [0.4; 0];
% X1_infeas = [0.4 0.4 0.4 0.15 ; 0 0.25 -0.25 0];

% if MULTIPLE_CIRCLES
%     X1_infeas = [X1_infeas, [-1 -1; -1 1]];
% X1_infeas = [X1_infeas, [-1; 1]];
%     
% end


% X1_infeas = [-1; -1];

opt.X = X;
if FEAS
    opt.X0 = X0_feas;
    opt.X1 = X1_feas;
else
    opt.X0 = X0_infeas;
    opt.X1 = X1_infeas;
end
% opt.verbose = 1;
order_range = [1, 4];
%% set up the problem
if SOLVE
order = 4;
d = 2*order;
% gamma = sdpvar(1,1);
[v, cv, mv] = polynomial([t;x],d);
[z1p, cz1p, mz1p] = polynomial([t;x],d);
[z2p, cz2p, mz2p] = polynomial([t;x],d);
[z1n, cz1n, mz1n] = polynomial([t;x],d);
[z2n, cz2n, mz2n] = polynomial([t;x],d);

v0 = replace(v, [t; x], [0; X0_infeas]);
% v1 = replace(v, [t;x], [1, X1_infeas]);
objective = -sum((x-X1_infeas).^2);
 
Lv = jacobian(v, t) + z1p + z2p + z1n + z2n;

%lie constraint
[pLie_m, consLie_m, GramLie_m] = psatz(Lv, X_moon, order, [t;x]);
[pLie_c, consLie_c, GramLie_c] = psatz(Lv, X_circ, order, [t;x]);

%objective constraint
[pobj_m, consobj_m, Gramobj_m] = psatz(v-objective, X_moon, order, [t;x]);
[pobj_c, consobj_c, Gramobj_c] = psatz(v-objective, X_circ, order, [t;x]);

%coefficients in zeta
coeff_1 = coefficients(z1p - z1n - jacobian(v, x(1)), [t;x]);
coeff_2 = coefficients(z2p - z2n - jacobian(v, x(2)), [t;x]);
con_coeff = [coeff_1==0; coeff_2==0];

%zeta nonnegativity

[pz1p_m, consz1p_m, Gramz1p_m] = psatz(z1p, X_moon, order, [t;x]);
[pz1p_c, consz1p_c, Gramz1p_c] = psatz(z1p, X_circ, order, [t;x]);

[pz2p_m, consz2p_m, Gramz2p_m] = psatz(z2p, X_moon, order, [t;x]);
[pz2p_c, consz2p_c, Gramz2p_c] = psatz(z2p, X_circ, order, [t;x]);

[pz1n_m, consz1n_m, Gramz1n_m] = psatz(z1n, X_moon, order, [t;x]);
[pz1n_c, consz1n_c, Gramz1n_c] = psatz(z1n, X_circ, order, [t;x]);

[pz2n_m, consz2n_m, Gramz2n_m] = psatz(z2n, X_moon, order, [t;x]);
[pz2n_c, consz2n_c, Gramz2n_c] = psatz(z2n, X_circ, order, [t;x]);

cons = [consLie_m; consLie_c; consobj_m; consobj_c; con_coeff;...
    consz1p_m; consz1p_c; ...
    consz2p_m; consz2p_c;...
    consz1n_m; consz1n_c;...
    consz2n_m; consz2n_c];

opts = sdpsettings('solver', 'mosek');

sol = optimize(cons, -v0, opts);

cost_rec = value(-v0);

% IM = set_manager(opt);
% order = 5;
% d = 2*order;
% out = IM.check_connected(d);
% out = IM.climb_connected(order_range);
end

%% plot the result
if SAMPLE
    
    [test, X_moon_func]=constraint_eval(X_moon, opt.x, [0;0]);
    [test, X_circ_func]=constraint_eval(X_circ, opt.x, [0;0]);
    
    
end
    

if PLOT
    bplot = moon_plotter(opt, out);
    bplot.inner_x = inner_x;
    bplot.inner_rad = inner_rad;
    bplot.circle = struct('center', circle_center, 'rad', circle_rad);
    bplot.set_plot();
%     bplot.contour_2d();
end