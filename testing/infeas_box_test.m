%try to do proofs of disconnectedness of double-lobe region

%This is a farkas certificate of infeasibility
%there exists a function v where v <= 0 on X0, v >= 0 on X1, and v
%decreases along trajectories. Trajectories will never cross the surface
%v(t,x) = 0, and this level sets separates X0 and X1.

    opt = set_path_options;
    opt.x = sdpvar(2,1);
    opt.Tmax = 2;

    opt.scale = 0;

%     order = 2;
%     order = 3;
    order = 4;
%     order = 5;
    d =2*order;
    T = 3; %maximum time

    X0_feas = [1.25; -1];
    X1_feas = [1.5; 0.5];


%     X0_infeas = [-0.75; 0];
            X0_infeas = [-0.75; 0.5];
%     X1_infeas = [1.5; 0.5];
%         X1_infeas = [1; 0.5];
%             X1_infeas = [1; -0.5];
            X1_infeas = [0.5; 0.75]; %needs order 4

    opt.box = 0;

    %constraint set
    f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);
    X.ineq = f(opt.x);
    X = fill_constraint(X);

    
    opt.time_indep = 0;
    
    opt.X = X;
    
    %matlab doesn't have deep copies.
    %annoying
    
    
    BOX = 1;
        
    %feasible test
    opt.X0 = X0_feas;
    opt.X1 = X1_feas;
    if BOX
        out_feas   = set_path_infeas_box(opt, order);        
    else
        out_feas   = set_path_infeas(opt, order);        
    end 
    %infeasible test    
    opt.X0 = X0_infeas;
    opt.X1 = X1_infeas;
    if BOX
        out_infeas   = set_path_infeas_box(opt, order);        
    else          
        out_infeas   = set_path_infeas(opt, order);        
    end

    disp(['Farkas Certificate: feas ', num2str(out_feas.farkas), ...
        ', infeas ', num2str(out_infeas.farkas)]);
    
    %should be feas 0, infeas 1
  