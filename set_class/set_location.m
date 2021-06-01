classdef set_location < handle
    %SET_LOCATION A cell of a partition (location) for set connectedness problems
    
    properties
        options = [];
        time_range = []; %range Tmin, Tmax of validity
        
        box = [];    %box in space
        poly = struct('v', [], 'zeta', [], 'nonneg', []);        
        
        id = 0;
        %other locations
        prev = {};
        next = {};
 
        %TODO: it is assumed for this implementation that X0 and X1 are
        %numeric. This should change later.
        
    end
    
    methods
        function obj = set_location(options, time_range, box, id)
            %SET_LOCATION Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
        
            obj.options = options;
            obj.time_range = time_range;
            obj.box = box;
            
            if nargin >= 4
                obj.id = id;
            end
            
            %adjacent locations
            n = length(options.x);            
            obj.prev = cell(n+1,1 );         
            obj.next = cell(n+1,1 );
           
               
          
        end
        
        function [prog_loc]= make_program(obj, d)
            
            %create program based purely on functions inside this location
            %ignore adjacency for now
            [poly_var, coeff_var] = obj.make_poly(d);
            
           
            
            [cons_loc, coeff_loc, nonneg_loc] = obj.make_cons(d, poly_var);
%             [cc_feas] = obj.make_cons(d, nonneg_feas);
            
            
%             cc_var_infeas = coef_con(coeff_var, []);
            
            prog_loc = struct('nonneg', nonneg_loc, 'poly', poly_var, ...
                'objective', 0, 'cons', cons_loc, 'coeff', [coeff_var; coeff_loc]);
        end
        
        function [poly_out, coeff_out] = make_poly(obj,d)
            %MAKE_POLY Create polynomials v and zeta for disconnectedness
            %certificate
            
            t = obj.options.t;
            x = obj.options.x;
            [v, cv] = polynomial([t; x], d);

            n = length(x);
            
            zeta = [];
            coeff_zeta = [];
            for i = 1:n
                [pzeta, czeta] = polynomial([t; x], d);
                zeta = [zeta; pzeta];
                coeff_zeta = [coeff_zeta; czeta];
            end
            
            %polynomials and derivatives for spline constraints
            v_vec = [v; jacobian(v, [t;x])'];            
            poly_vec = [v_vec; zeta];
            
            poly_out=struct('v', v, 'zeta', zeta, 't', t, 'x',x, 'vec', poly_vec);
            coeff_out = [cv; coeff_zeta];
        end

        
        %% Generate constraints within location
        function [cons_loc, coeff_loc, nonneg_loc, poly_var] = make_cons(obj, d, poly_var)
            %generate constraints within this location
            
            if nargin < 3
                poly_var = obj.make_poly(d);
            end
            
            [con_init, nonneg_init] = obj.make_init_con(poly_var);
            [con_term, nonneg_term] = obj.make_term_con(poly_var);

            [nonneg_loc] = obj.make_lie_nonneg(poly_var);
           
            [con_lie, coeff_lie] = obj.make_lie_con(d, nonneg_loc);
             
            %package up constraints and nonnegative functions
            cons_loc = [con_init; con_term; con_lie];
            coeff_loc = coeff_lie;
            
            nonneg_loc.init = nonneg_init;
            nonneg_loc.term = nonneg_term;
            
        end
        
        function [con_init, nonneg_init] = make_init_con(obj, poly_var)
            %constraints on initial points X0
            con_init = [];
            nonneg_init = [];

            %only form initial constraints if the valid time range
            %includes the beginning time

            if obj.time_range(1) == 0
                X0 = obj.options.X0;
                X0_contained = all((X0 >= obj.box(:, 1)) .* (X0 <= obj.box(:, 2)),1);

                X0_active = find(X0_contained);

                con_init = [];
                if isempty(X0_active)                
                    nonneg_init = [];
                else
                    t = poly_var.t;
                    x = poly_var.x;  
%                     v = poly_var.v;
                    v0 = replace(poly_var.v, t, 0);
                    nonneg_init = v0 - 1; %v0 >= 1

                    for i = 1:length(X0_active)
                        v0_val = replace(v0, x, X0(:, i));
                        con_init = [con_init; v0_val >= 0];

                    end
                    con_init = con_init:['Cell ', num2str(obj.id), ' initial'];
                end
                
                
            end
        end
   
        function [con_term, nonneg_term] = make_term_con(obj, poly_var)
                %constraints on terminal points X1
            
                con_term = [];
                nonneg_term = [];
                
                %only form terminal constraints if the valid time range
                %includes the end time
                
                if obj.time_range(2) == obj.options.Tmax
                    X1 = obj.options.X1;
                    X1_contained = all((X1 >= obj.box(:, 1)) .* (X1 <= obj.box(:, 2)),1);

                    X1_active = find(X1_contained);


                    if ~isempty(X1_active)                
                        t = poly_var.t;
                        x = poly_var.x;  
                        v = poly_var.v;
                        if obj.options.scale
                            v1 = replace(v, t, 1);
                        else
                            v1 = replace(v, t, obj.options.Tmax);
                        end
                        nonneg_term = -v1; %v1 <= 1
                        for i = 1:length(X1_active)
                            v1_val = replace(v1, x, X1(:, i));
                            con_term = [con_term; v1_val >= 0];
                        end
                        con_term = con_term:['Cell ', num2str(obj.id), ' terminal'];
                    end
                end
            end
        
  
        function [nonneg_lie] = make_lie_nonneg(obj, poly_var)
            %constraints enforcing that v increases along controlled
            %trajectories
            
            
            %time-scaling
            if obj.options.scale
                scale_weight = obj.options.Tmax;
            else
                scale_weight = 1;
            end
            
            
            t = obj.options.t;
            x = obj.options.x;
            n = length(x);
            
            v = poly_var.v;
            zeta = poly_var.zeta;
            zeta_sum = sum(zeta);
            
            %occupation measures
            Lv = jacobian(v, t) - scale_weight*jacobian(v,x)*ones(n, 1) - zeta_sum;
            nonneg_lie.occ = Lv; %Lv >= 0
            
            %box-input measures
            term_ui = [];
            for i = 1:n
                term_ui = [term_ui; zeta(i) + scale_weight*2*jacobian(v,x(i))];                
            end            
            
            nonneg_lie.u = term_ui;     
            
            nonneg_lie.slack = zeta;
            
        end
        
        function [con_lie, coeff_lie] = make_lie_con(obj, d, nonneg)

            %lie constraints and control/slack constraints
            
            %variables           
            t = obj.options.t;
            x = obj.options.x;
            n = length(x);  
            
             X_cell  = obj.prep_space_cell(obj.options.X);            
            
            
                        %occupation measure (+box, complement)
            if obj.options.scale
                Tsupp = t*(1-t);
            else
                Tsupp = t*(obj.options.Tmax - t);
            end   
            
            con_occ = [];
            coeff_occ = [];
            
            con_u = []; 
            coeff_u = []; 
            
            con_slack = [];
            coeff_slack = [];
            
            for i = 1:length(X_cell)
                X_curr = X_cell{i};
                X_curr.ineq = [Tsupp; X_curr.ineq];
                
                %occupation
                [con_occ_curr, coeff_occ_curr] = obj.make_psatz(d, X_curr, nonneg.occ, [t; x]);
                con_occ_curr = con_occ_curr:['Cell ', num2str(obj.id), ' Lie base '];
                con_occ = [con_occ; con_occ_curr];
                coeff_occ = [coeff_occ; coeff_occ_curr];
                
                
                for j = 1:n
                    %box 
                    [con_u_curr, coeff_u_curr] = obj.make_psatz(d, X_curr, nonneg.u(j), [t; x]);
                    con_u_curr = con_u_curr:['Cell ', num2str(obj.id), ' Lie input ', num2str(j)];
                    con_u = [con_u; con_u_curr];
                    coeff_u = [coeff_u; coeff_u_curr];

                
                    
                    %complement
                    [con_slack_curr, coeff_slack_curr] = obj.make_psatz(d, X_curr, nonneg.slack(j), [t; x]);
                    con_slack_curr = con_slack_curr:['Cell ', num2str(obj.id), ' Slack input ', num2str(j)];
                    con_slack = [con_slack; con_slack_curr];
                    coeff_slack = [coeff_slack; coeff_slack_curr];
                
                end
                
                
            end    
            
            con_lie = [con_occ; con_u; con_slack];
            coeff_lie = [coeff_occ; coeff_u; coeff_slack];
            
        end
                %% adjacency constraints
                
%          %this may be appropriate to be implemented by the manager
%         function [con_adj] = make_adjacency_con(obj)
%             %create equality constraints between polynomials on boundaries
%             %of neighboring cells
%             con_adj = [];
%             t = obj.poly.t;
%             x = obj.poly.x;
%             vars = [t; x];
%             
%             for i = 1:length(obj.next)
%                 if ~isempty(obj.next{i})  
%                     %get boundary polynomials
%                     pv_curr = obj.get_adjacency_poly(i-1, 1);
%                     pv_next = obj.next{i}.get_adjacency_poly(i-1, 0);
%                     
%                     %set boundary polynomials equal to each other
%                     coeff_pv = coefficients(pv_curr - pv_next, vars);
%                     con_adj = [con_adj; coeff_pv == 0];
%                 end
%             end
%         end
        

        
        function poly_vec = get_adjacency_poly(obj,poly_var, ind, direction )
            %evaluate polynomials on an edge of the current box
            %used for the spline constraints over the partitions
            
            %Inputs
            %   poly_var:   container for polynomials
            %   ind:        index of boundary dimension (0: time, 1+: space)
            %   direction:  previous = 0 (default), next = 1
            
            
            if nargin < 4
                direction = 0;
            end
            t = poly_var.t;
            x = poly_var.x;

            vec = poly_var.vec;
            if ind == 0
                t_val = obj.time_range(1+direction);
                
                if obj.options.scale
                    t_val = t_val/obj.options.Tmax;
                end
                poly_vec = replace(vec, t, t_val);
            else
                x_val = obj.box(ind, 1+direction);
                poly_vec = replace(vec, x(ind), x_val);
            end
            
        end
            

        
        
        function [cons, coeff] = make_psatz(obj, d, X, f, vars)
            %MAKE_PSATZ a positivestellensatz expression to impose that
            %f(vars) >= 0 on the region X
            if isstruct(X)
                %X1 is a set, nonnegativity of region
%                 set1 = true;
                X = fill_constraint(X);

                [p, cons, coeff] = constraint_psatz(f, X, vars, d);        
            else  
                %X is a point, nonnegativity of evaluation
                f_pt = replace(f, vars, X);
                cons = f_pt>= 0; 
                coeff = [];
            end
            
%             cc_curr = coef_con(coeff, cons);
        end
        
            function X_cell = prep_space_cell(obj, X)
            %PREP_SPACE_CELL: Split up X into cells
            if iscell(X)
                %X is already a cell, nothing to be done
                X_cell = cell(length(X), 1);
                for i = 1:length(X)
                    X_cell{i} = fill_constraint(X{i});
                end
            else
                if isstruct(X)
                    %wrap up the structure in a cell
                    X_cell = {X};
                else
                    %X is supported at discrete points
                    %Each column of X is a possible origin datapoint
                    if size(X, 2) == 1
                        X_cell = {X};
                    else
                        Xt = X';
                        Xt_cell = mat2cell(Xt, ones(size(X, 2), 1));
                        X_cell = cellfun(@(x) x', Xt_cell, 'UniformOutput', false);
                    end
                end
            end
        end
               
     
        %% overloads
  
        
    end
end

