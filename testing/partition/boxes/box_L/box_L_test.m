SOLVE =1;
PLOT = 0;
SAMPLE = 0;
EVAL = 0;

rng(55, 'twister');

%ring-circle test exhibits circular symmetry, depends only on the radius
%solve a 1d set connectedness program based on the radius 

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
% opt.Tmax = 1;
opt.box = [0, 1; 0, 1];

opt.scale = 1;
opt.verbose = 1;



%% support set (in unit ball)


%two L shapes
boxes = { [0.2, 0.3;0.3, 0.7], ...
          [0.2, 0.6;0.6, 0.7], ...
          [0.4, 0.8;0.3, 0.4], ...
          [0.7, 0.8; 0.3, 0.7]};

% boxes = { [0.2, 0.3;0.3, 0.7], ...
%           [0.2, 0.6;0.6, 0.7], ...
%           [0.4, 0.8;0.3, 0.4]};

X = cell(length(boxes), 1);
boxfunc = cell(length(boxes), 1);
Tbox = zeros(length(boxes), 1);

n = length(opt.x);
%create constraints for each box
for i = 1:length(boxes)
    boxcurr = boxes{i};
    Xcurr = struct('ineq', [opt.x-boxcurr(:, 1); boxcurr(:, 2) - opt.x], 'eq', []);
    X{i} = Xcurr;
    
    %also numerical support evaluations
    boxfunc{i} = @(x) all((x - boxcurr(:, 1)) > 1e-8) && all((x - boxcurr(:, 2)) < -1e-8);
    
    Tbox(i) = norm(diff(boxcurr'), n);
end

Tmax = sum(Tbox);

X_func = @(x) any(cellfun(@(b) b(x), boxfunc));
      
      
%% initial sets
X0 = [0.25; 0.65];
X1 = [0.75; 0.35];

% POINTS0 = 1;\
% POINTS1 = 1;


% X0 = [0.8];
% X1 = [0.3];
% 
% X_in = struct('ineq', [R_circ - x; x] , 'eq', []);
% X_out = struct('ineq', [R_ring_outer-x; x-R_ring_inner], 'eq',[]);

% X = {X_in, X_out};

opt.X = X;
opt.X0 = X0;
opt.X1 = X1;
opt.epsilon = 0.01;

opt.Tmax = Tmax;

if SOLVE
% IM = set_manager(opt);

% spacing = [1; 1];

% spacing = [1; 2; 2];
spacing = [4; 1; 1];

% spacing = [2; 2; 2];
% spacing = [4;2];

% spacing = [5; 6];
% spacing= [3; 3; 3];

% spacing = [10; 4; 4];
SM = set_manager_partition(opt, spacing);


% order = 5;
order = 4;
% order = 3;
d = 2*order;

out = SM.check_connected(d);
% out.farkas
end


if SAMPLE
    
    supp_func = @(t,x) support_event(t, x, X_func);

    s_opt = set_sample_options;
    % x0 = @() [0;0];

    s_opt.x0 = @() X0;

    s_opt.Tmax = opt.Tmax;
    s_opt.dt = 0.1;
%     s_opt.dt = 0.025;
    s_opt.X_func = supp_func;

    Np = 100;
%     Np = 300;

    out_sim=set_walk(Np, s_opt);

end

if (SAMPLE || EVAL) && (out.status == conn_status.Disconnected)
  
    out_sim = set_traj_eval(out, out_sim);

end

% 
% %TODO: write plotting code for ellipses
if PLOT &&  out.status == conn_status.Disconnected
%     bplot = ring_circ_plotter(opt, out);
    bplot = box_plotter(opt, out, out_sim);
    bplot.boxes = boxes;

    bplot.traj_2d();
    bplot.v_plot();
    bplot.nonneg_traj();
    bplot.nonneg_zeta();
        
%     bplot.contour_2d();
end
