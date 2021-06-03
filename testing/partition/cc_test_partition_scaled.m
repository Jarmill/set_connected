SOLVE = 1;
SAMPLE = 1;
EVAL = 1;
PLOT = 1;

%this example works without partitioning

opt = set_path_options;

opt.t = sdpvar(1, 1);
opt.x = sdpvar(2,1);
% opt.Tmax = 20;
% opt.Tmax = 10;
opt.Tmax = 5;
% opt.Tmax = 2;
% opt.Tmax = 151;
opt.verbose = 1;

opt.epsilon = 0.01;

opt.scale = 1;

fscale =  @(varnew)(-0.749969140625+1.4666015625.*varnew(1, :)+1.4761.*varnew(2, :)-12.3596191406.*varnew(1, :).^4-4.42050625.*varnew(2, :).^4+10.2172851562.*varnew(1, :).^2-0.2417875.*varnew(2, :).^2+3.9421875.*varnew(1, :).*varnew(2, :).^2+1.0875.*varnew(1, :).*varnew(2, :)-3.2958984375.*varnew(1, :).^3-2.4389.*varnew(2, :).^3);

% figure
fimplicit(@(x,y)fscale([x;y]), [-1 1 -1 1]);
X=struct;
X.ineq = fscale(opt.x);
X = fill_constraint(X);


opt.X = X;



opt.X0 = [-0.466666666666667,-0.466666666666667,-0.466666666666667,-0.573333333333333,-0.573333333333333,-0.573333333333333;0.0689655172413793,-0.0689655172413793,-0.206896551724138,0.0689655172413793,-0.0689655172413793,-0.206896551724138];
opt.X1 = [0.733333333333333,0.733333333333333,0.733333333333333,0.626666666666667,0.626666666666667,0.626666666666667;0.206896551724138,0.0689655172413793,-0.0689655172413793,0.206896551724138,0.0689655172413793,-0.0689655172413793];

order_range = [1, 3];

opt.box = 1;

if SOLVE
% IM = set_manager(opt);
% out = IM.check_connected(d);

spacing= [1; 2; 1];
% spacing = [2;2;1];
% spacing = [4; 2; 1];

% spacing = [10; 4; 4];
SM = set_manager_partition(opt, spacing);

out = SM.check_connected(6);
% out = SM.climb_connected(order_range);
end


if SAMPLE
    %     Np = 100;
        [test, X_func]=constraint_eval(X, opt.x, [0;0]);

        supp_func = @(t,x) support_event(t, x, X_func);

    s_opt = set_sample_options;
    % x0 = @() [0;0];

    s_opt.x0 = @() opt.X0(:,randi( size(opt.X0, 2)));

    
    s_opt.u_boundary = 0;
    s_opt.Tmax = opt.Tmax;
    s_opt.dt = 0.05;
%     s_opt.dt = 0.1;
    s_opt.X_func = supp_func;
%     s_opt.nonneg_func = out.func.nonneg;

    Np = 30;

    out_sim=set_walk(Np, s_opt);

end


if SAMPLE || EVAL
   
    
    out_sim = set_traj_eval(out, out_sim);
end



if PLOT && out.status == conn_status.Disconnected
    lobe_plot = lobe_plotter(opt, out, out_sim);

    lobe_plot.v_plot();
    lobe_plot.nonneg_traj();
    lobe_plot.nonneg_zeta();
    
%         lobe_plot.contour_2d();
    
%     lobe_plot.contour_3d();
end