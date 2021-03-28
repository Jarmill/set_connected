classdef infeas_manager
    %INFEAS_MANAGER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        options = [];
        poly = struct('v', [], 'zeta', []);        
        
        cc = coef_con([], []);

    end
    
    methods
        function obj = infeas_manager(options)
            %INFEAS_MANAGER Construct an instance of this class
            %   Detailed explanation goes here
            obj.options = options;
            if isempty(obj.options.t)
                obj.options.t = sdpvar(1, 1);
            end
        end
        
        function [cc_all, poly_var]= make_program(obj, d)
            [poly_var, coeff_var] = obj.make_poly(d);
            
            nonneg = obj.form_nonneg(poly_var);
            
            
            [cc_cons] = obj.make_cons(d, nonneg);
            
            cc_var = coef_con(coeff_var, []);
            
            cc_all = [cc_var; cc_cons];
        end
        
        %% variables
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
            
            poly_out=struct('v', v, 'zeta', zeta, 't', t, 'x', x);
            coeff_out = [cv; coeff_zeta];
        end
        
        %% Constraints
        
        function nonneg = form_nonneg(obj, poly)
            %FORM_NONNEG create functions that should be nonnegative
            t = obj.options.t;
            x = obj.options.x;
            n = length(x);
            
            %polynomials
            v = poly.v;
            zeta = poly.zeta;
            zeta_sum = sum(zeta);
            
            %time-scaling
            if obj.options.scale
                scale_weight = obj.options.Tmax;
            else
                scale_weight = 1;
            end
            
            nonneg = struct('init', [], 'term', [], 'occ', [], 'u', [], 'slack', []);
            %init:      initial measure
            %term:      terminal measure
            %occ:       occupation measure
            %u:         box-input occupation measure
            %slack:     box-input complement occupation measure (abs. cont.)            
            
            
            %TODO: make sensible choices of signs
            
            %initial set
            v0 = replace(v, t, 0);
            nonneg.init = v0 - 1; %v0 >= 1
            
            %terminal set
            if obj.options.scale
                vT = replace(v, t, 1);
            else
                vT = replace(v, t, obj.options.Tmax);
            end
            nonneg.term = -vT - 1; %vT <= -1
            
            %occupation measures
            Lv = jacobian(v, t) - scale_weight*jacobian(v,x)*ones(n, 1) - zeta_sum;
            nonneg.occ = Lv; %Lv >= 0
            
            %box-input measures
            term_ui = [];
            for i = 1:n
                term_ui = [term_ui; zeta(i) + scale_weight*2*jacobian(v,x(i))];                
            end            
            
            nonneg.u = term_ui;     
            
            nonneg.slack = zeta;
        end             
        
        function cc_curr = make_psatz(obj, d, X, f, vars)
            %f >= 0 on X
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
            
            cc_curr = coef_con(coeff, cons);
        end
        
        function cc_all = make_cons(obj, d, nonneg)
            %make constraints for set disconnectedness program
            
            %sets
            X0_cell = obj.prep_space_cell(obj.options.X0);
            X1_cell = obj.prep_space_cell(obj.options.X1);
            X_cell  = obj.prep_space_cell(obj.options.X);            
            
            %variables           
            t = obj.options.t;
            x = obj.options.x;
            n = length(x);                                    
           
            %initial set
            cc_init = coef_con();
            for i = 1:length(X0_cell)
                cc_init_curr = obj.make_psatz(d, X0_cell{i}, nonneg.init, x);
                cc_init = cc_init.vertcat(cc_init_curr);
            end
             
            %terminal set
            cc_term = coef_con();
            for i = 1:length(X1_cell)
                cc_term_curr = obj.make_psatz(d, X1_cell{i}, nonneg.term, x);
                cc_term = [cc_term; cc_term_curr];
            end
             

            %occupation measure (+box, complement)
            if obj.options.scale
                Tsupp = t*(1-t);
            else
                Tsupp = t*(obj.options.Tmax - t);
            end            
            
            cc_occ = coef_con();
            cc_u = coef_con();
            cc_slack = coef_con();
            for i = 1:length(X_cell)
                X_curr = X_cell{i};
                X_curr.ineq = [Tsupp; X_curr.ineq];
                
                %occupation
                cc_occ_curr = obj.make_psatz(d, X_curr, nonneg.occ, [t; x]);
                cc_occ = [cc_occ; cc_occ_curr];
                
                
                for j = 1:n
                    %box 
                    cc_u_curr = obj.make_psatz(d, X_curr, nonneg.u(j), [t; x]);
                    cc_u = [cc_u; cc_u_curr];
                    
                    %complement
                    cc_slack_curr = obj.make_psatz(d, X_curr, nonneg.slack(j), [t; x]);
                    cc_slack = [cc_slack; cc_slack_curr];
                end
                
                
            end    
            
            %pack everything up
%             cc_all = [cc_init; cc_term; cc_occ; cc_u; cc_slack];
            cc_all = [cc_init; cc_term];
            cc_all = [cc_all; cc_occ];
            cc_all = [cc_all; cc_u];
            cc_all = [cc_all; cc_slack];
            
        end
        
        
        %% constraint helper functions
        
        function X_cell = prep_space_cell(obj, X)
            %PREP_SPACE_CELL: Split up X into cells
            if iscell(X)
                %X is already a cell, nothing to be done
                X_cell = X;
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
        
        
        function cc_init = con_init(obj, X0, d)
%             X0 = obj.options.X0;
            t = obj.options.t;
            x = obj.options.x;
            n = length(x);
            v = obj.poly.v;
            
            %value function at (0, X0) is positive
            if isstruct(X0)    
                %X0 is a set                
                X0 = fill_constraint(X0);
                v0 = replace(v, t, 0);
                [p0, cons0, coeff0] = constraint_psatz(v0 - 1, X0, x, d);
            else
%                 X0 is a point                    
                v0 = replace(v, [t;x], [0; X0]);
%                 cons0 = (v0 == 1);
                cons0 = (v0 -1 >= 0);
                coeff0 = [];    
            end
            
            cc_init = coef_con(coeff0, cons0);                                    
        end
        
        function cc_term = con_term(obj, X1, d)
            %value function on (T, X1) is negative
            t = obj.options.t;
            x = obj.options.x;
            n = length(x);
            v = obj.poly.v;
            
            if isstruct(X1)
                %X1 is a set
%                 set1 = true;
                X1 = fill_constraint(X1);
                if obj.options.scale
                    vT = replace(v, t, 1);
                else
                    vT = replace(v, t, obj.options.Tmax);
                end

                [pT, consT, coeffT] = constraint_psatz(-vT -1, X1, x, d);        
            else
    %             X1 is a point                
                if obj.options.scale
                    vT = replace(v, [t; x], [1; X1]);
                else
                    vT = replace(v, [t; x], [obj.options.Tmax; X1]);
                end
%                 [pT, consT, coeffT] = constraint_psatz(vT, X1, x, d);
%                 consT = (vT == -1);
                consT = -vT -1>= 0; 
                coeffT = [];
            end
            
            cc_term = coef_con(coeffT, consT);
        end
        
        function [cc_occ, cc_u] = con_occ(obj, X_region, d)
            %value function v increases along all controlled trajectories
            %therefore impossible for trajectories to go from X0 to X1
            %variables
            t = obj.options.t;
            x = obj.options.x;
            n = length(x);
            
            %polynomials
            v = obj.poly.v;
            zeta = obj.poly.zeta;
            zeta_sum = obj.poly.zeta_sum;
            
            if obj.options.scale
                scale_weight = obj.options.Tmax;
            else
                scale_weight = 1;
            end
            
            %occupation measure

            Lv = jacobian(v, t) - scale_weight*jacobian(v,x)*ones(n, 1) - zeta_sum;

            [pL, consL, coeffL] = constraint_psatz(Lv, X_region, [t; x], d);
            
            cc_occ = coef_con(coeffL, consL);
            
            %box measure
            cons_u = [];
            coeff_u = [];
            for i = 1:n
                term_ui = zeta(i) + scale_weight*2*jacobian(v,x(i));
                [p_ui, cons_ui, coeff_ui] = constraint_psatz(term_ui, X_region, [t; x], d);

                cons_u = [cons_u; cons_ui];
                coeff_u = [coeff_u; coeff_ui];
            end
            
            cc_u = coef_con(coeff_u, cons_u);
        end
        
        
        
        
    end
end

