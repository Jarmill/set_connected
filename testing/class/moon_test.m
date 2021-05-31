%moon test 
%there was a bug in the code, X was set to {} before adding circles
%now the certificate fails. Will try partitioning


SOLVE = 1;
PLOT = 1;
SAMPLE = 0;
EVAL = 0;
FEAS = 0;

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
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
    X_circ.ineq = rx^2 - (opt.x(1)-cx(1))^2 - (opt.x(2)-cx(2))^2;
    X = [X; {X_circ}];
end


X0_feas = [0; 0.9];
X1_feas = [0; -0.9];

% X0_infeas = [0 0 0.7 0.7 -0.8; 0.8 -0.8 0.65 -0.65 0 ];
X0_infeas = X0_feas;
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

if SOLVE
IM = set_manager(opt);
order = 5;
d = 2*order;
out = IM.check_connected(d);
% out = IM.climb_connected(order_range);
end

if PLOT
    bplot = moon_plotter(opt, out);
    bplot.inner_x = inner_x;
    bplot.inner_rad = inner_rad;
    bplot.circle = struct('center', circle_center, 'rad', circle_rad);
%     bplot.set_plot();
    bplot.contour_2d();
end