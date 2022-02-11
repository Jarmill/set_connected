%moon test 
%there was a bug in the code, X was set to {} before adding circles
%now the certificate fails. Will try partitioning


SOLVE = 1;
PLOT = 1;
SAMPLE = 0;
EVAL = 0;
FEAS = 0;

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('b', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;
vars.b = b;
lsupp = loc_support(vars);
lsupp = lsupp.set_box(1);
% opt = set_path_options;

% opt.t = sdpvar(1, 1);
% vars.x = sdpvar(2,1);
% opt.Tmax = 2;
lsupp.Tmax = 2;

%controlled dynamics
f = 2*b-1;

%moon set
inner_rad = 0.7;
inner_x = 0.35;

% X_moon = struct;
X_moon= [1 - vars.x(1)^2 - vars.x(2)^2; 
                (vars.x(1)-inner_x)^2 + vars.x(2)^2 - inner_rad^2];
% X_moon = fill_constraint(X_moon);

X = {X_moon>=0};

%circles
MULTIPLE_CIRCLES = 0;

if MULTIPLE_CIRCLES
    circle_rad = [0.4; 0.3; 0.3];
    circle_center = [0.4, 0;
                 -1, 1;
                 -1, -1]';

%     circle_rad = [0.4; 0.3];
% circle_center = [0.4, 0;
%                  -1, 1]';

else    
    circle_rad = [0.4];
    circle_center = [0.4, 0]';
end




% circle_rad = [0.4; 0.3; 0.3];
% circle_center = [0.4, 0;
%                  -1, 1;
%                  -1, -1]';

% circle_rad = [1; 1]*0.3;
% circle_center = [-1, 1;
%                  -1, -1]';

X_circ = {};
for i = 1:length(circle_rad)
%     X_circ = struct;
    rx = circle_rad(i);
    cx = circle_center(:, i);
    X_circ= rx^2 - (vars.x(1)-cx(1))^2 - (vars.x(2)-cx(2))^2;
    X = [X; {X_circ>=0}];
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

lsupp.X = X;
if FEAS
    lsupp.X_init = X0_feas;
    objective = -sum((x-X1_feas).^2);
else
    lsupp.X_init = X0_infeas;
    objective = -sum((x-X1_infeas).^2);
end
% opt.verbose = 1;
order_range = [1, 4];

if SOLVE
% IM = set_manager(opt);
PM = peak_manager(lsupp, f, objective);

order = 3;
% order = 5;
d = 2*order;
% out = IM.check_connected(d);
sol = PM.run(order);
disp(sol.obj_rec)
% out = IM.climb_connected(order_range);
end

% if SAMPLE
%     
%     [test, X_moon_func]=constraint_eval(X_moon, vars.x, [0;0]);
%     [test, X_circ_func]=constraint_eval(X_circ, vars.x, [0;0]);
%     
%     
% end
    

% if PLOT
%     bplot = moon_plotter(opt, out);
%     bplot.inner_x = inner_x;
%     bplot.inner_rad = inner_rad;
%     bplot.circle = struct('center', circle_center, 'rad', circle_rad);
%     bplot.set_plot();
% %     bplot.contour_2d();
% end