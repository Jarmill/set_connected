classdef infeas_manager
    %INFEAS_MANAGER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        options = [];
        poly = struct('v', [], 'zeta', [], 'zeta_sum', []);        
        
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
        
        function [obj] = make_program(obj, d)
            [obj.poly, cc_poly] = obj.make_poly(d);
            [cc_cons] = obj.make_cons(d);
            
            
            obj.cc = cc_poly.append(cc_cons);
            %solve yalmip
        end
        
        function [poly_out, cc_out] = make_poly(obj,d)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
            t = obj.options.t;
            x = obj.options.x;
            [v, cv] = polynomial([t; x], d);

            n = length(x);
            
            zeta = [];
            coeff_zeta = [];
            cons_zeta_sos = [];
            zeta_sum = 0;
            for i = 1:n
                [pzeta, czeta] = polynomial([t; x], d);
                coeff_zeta = [coeff_zeta; czeta];
                cons_zeta_sos = [cons_zeta_sos; sos(pzeta)];

                zeta = [zeta; pzeta];
                zeta_sum = zeta_sum + pzeta;
            end
            
            poly_out=struct('v', v, 'zeta', zeta, 'zeta_sum', zeta_sum, 't', t, 'x', x);
            coeff_out = [cv; coeff_zeta];
            cons_out = cons_zeta_sos;
            cc_out = coef_con(coeff_out, cons_out);
        end
        
        function nonneg = form_nonneg(obj, poly)
            %FORM_NONNEG %create functions that should be nonnegative
            
            
            nonneg = struct('init', [], 'term', [], 'occ', [], 'u', []);
            
        end
        
        function cc_all = make_cons(obj, d)
            %make constraints for set disconnectedness program
            
            X0_cell = obj.prep_space_cell(obj.options.X0);
            X1_cell = obj.prep_space_cell(obj.options.X1);
            X_cell  = obj.prep_space_cell(obj.options.X);            
            
            cc_all = coef_con([], []);
            
            for i = 1:length(X0_cell)
                cc_init_curr = obj.con_init(X0_cell{i}, d);
                cc_all = cc_all.append(cc_init_curr);
            end
            
            for i = 1:length(X1_cell)
                cc_term_curr = obj.con_term(X1_cell{i}, d);
                cc_all = cc_all.append(cc_term_curr);
            end
            
            for i = 1:length(X_cell)
                [cc_occ_curr, cc_u_curr] = obj.con_occ(X_cell{i}, d);
                cc_all = cc_all.append(cc_occ_curr);
                cc_all = cc_all.append(cc_u_curr);
            end            
            
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
                        Xt_cell = mat2cell(Xt', ones(size(Xt, 2), 1));
                        X_cell = cellfun(@(x) x', Xt_cell);
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

