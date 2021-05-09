classdef set_manager
    %SET_MANAGER 
    %Given sets X0 and X1 within a base space X,
    %determine if there exists a path between points in X0 and X1 (feas)
    %or there does not exist such a path (infeas)
    
    properties
        options = [];
        poly = struct('v', [], 'zeta', []);        
        
        cc = coef_con([], []);

    end
    
    methods
        function obj = set_manager(options)
            %SET_MANAGER Construct an instance of this class
            %   Detailed explanation goes here
            obj.options = options;
            if isempty(obj.options.t)
                obj.options.t = sdpvar(1, 1);
            end
        end
        
        %% program and recovery
        
        function out = climb_connected(obj, order_range)
            %meta-algorithm: increase degree between d_range(1):d_range(2)
            
            
            if nargin == 1
                order_range = [1; 5];
            end
            
            if length(order_range) == 1
                order_range = [1; d];
            end
            d_range = 2*order_range;
            
            for d_curr = d_range(1):d_range(2)
                out = obj.check_connected(d_curr);
                fprintf('order=%d, \t status=%s\n', d_curr/2, char(out.status))
                if (out.status ~= conn_status.Indeterminate)
                    break
                end
            end
            
        end
        
        function out = check_connected(obj, d)

            %the previously described 'prog_feas' feasibility program is an
            %outer approximation of the minimum time control, and is a 
            % misnomer. This program will return a time bound of Inf if the
            % two sets are disconnected. Other types of programs are
            % required to certify feasibility (density?)
            
            prog_infeas = obj.make_program(d);
            out_infeas = solve_program(obj, prog_infeas);
            farkas= ~out_infeas.problem;
            if farkas
%                 out_feas = struct('');
                status = conn_status.Disconnected;
            else
                status = conn_status.Indeterminate;          
                
%                 connected = 0;
%             else
%                 out_feas = solve_program(obj, prog_feas);
%                 connected = ~out_feas.problem;
            end
            
%             if farkas
%                 status = conn_status.Disconnected;
%             elseif connected
%                 status = conn_status.Connected;                
%             else
%                 status = conn_status.Indeterminate;
%             end
            
%             out = struct('status', status, 'infeas', out_infeas, 'feas', out_feas);
%             out = struct('status', status, 'infeas', out_infeas);
%             out = status;
            out = out_infeas;
            out.status = status;

            %else: it is inconclusive whether the sets are disconnected
            %out.farkas = 0 (default)
        end
        
%         function [prog_infeas, prog_feas]= make_program(obj, d)
        function [prog_infeas]= make_program(obj, d)
            [poly_var, coeff_var] = obj.make_poly(d);
            
            [nonneg_infeas, nonneg_feas] = obj.form_nonneg(poly_var);
            
            
            [cc_infeas] = obj.make_cons(d, nonneg_infeas);
%             [cc_feas] = obj.make_cons(d, nonneg_feas);
            
            
            cc_var_infeas = coef_con(coeff_var, []);
%             cc_var_feas = coef_con([poly_var.gamma; coeff_var], []);
            
            cc_all_infeas = [cc_var_infeas; cc_infeas];
%             cc_all_feas = [cc_var_feas; cc_feas];
            
%             
%             feas_objective = -poly_var.gamma;
%             feas_objective = 0;
%             prog_feas = struct('nonneg', nonneg_feas, 'poly', poly_var, ...
%                 'objective', feas_objective, 'cc', cc_all_feas);
            prog_infeas = struct('nonneg', nonneg_infeas, 'poly', poly_var, ...
                'objective', 0, 'cc', cc_all_infeas);
        end
        
        function [out] = solve_program(obj, prog)
            %solve SOS program in YALMIP, return solution    
            
            opts = sdpsettings('solver', obj.options.solver, 'verbose', obj.options.verbose);
            opts.sos.model = 2;
            
            [sol, monom, Gram, residual] = solvesos(prog.cc.con, prog.objective, opts, prog.cc.coef);
            
            out = struct('poly', [], 'problem', sol.problem, 'sol', [], 'block', [], 'func', []);
            if sol.problem == 0
                %the sets X0 and X1 are disconnected in time range [0, T]
                [out.poly, out.func] = obj.recover_poly(prog.poly, prog.nonneg);
                out.sol = sol;   
                out.block = struct;
                out.block.monom = monom;
                out.block.Gram = Gram;
                out.block.residual = residual;                            
            end
        end
        
        function [poly_eval, func_eval] = recover_poly(obj, poly_var, nonneg)
            %recover polynomials from computed disconnectedness certificate
            t = poly_var.t;
            x = poly_var.x;
            n = length(poly_var.x);
            
            %solved coefficients of v and zeta
            [cv,mv] = coefficients(poly_var.v,[poly_var.t; poly_var.x]);
            v_eval = value(cv)'*mv;

            [cz, mz] = coefficients(poly_var.zeta,[poly_var.t; poly_var.x]);
            if n == 1
                zeta_eval = value(cz)'*mz;
            else
                zeta_eval = value(cz)*mz;
            end
            
            %evaluations of v at initial and terminal times
            v0 = replace(v_eval, t, 0);            
            
            %terminal set
            if obj.options.scale
                v1 = replace(v_eval, t, 1);
            else
                v1 = replace(v_eval, t, obj.options.Tmax);
            end
                        
            poly_eval = struct('v', v_eval, 'zeta', zeta_eval, 'v0', v0, 'v1', v1);
            
            % form functions using helper function 'polyval_func'
            func_eval = struct;
            func_eval.v = polyval_func(v_eval, [t; x]);
            func_eval.zeta = polyval_func(zeta_eval, [t; x]);
            
            func_eval.v0 = polyval_func(v0, [x]);            
            func_eval.v1 = polyval_func(v1, [x]);
                                                            
            %TODO: function evaluations
        end
        
        
        %% variables
        function [poly_out, coeff_out] = make_poly(obj,d)
            %MAKE_POLY Create polynomials v and zeta for disconnectedness
            %certificate
            
            t = obj.options.t;
            x = obj.options.x;
            [v, cv] = polynomial([t; x], d);
            gamma = sdpvar(1,1);

            n = length(x);
            
            zeta = [];
            coeff_zeta = [];
            for i = 1:n
                [pzeta, czeta] = polynomial([t; x], d);
                zeta = [zeta; pzeta];
                coeff_zeta = [coeff_zeta; czeta];
            end
            
            poly_out=struct('v', v, 'zeta', zeta, 't', t, 'x', x, 'gamma', gamma);
            coeff_out = [cv; coeff_zeta];
        end
        
        %% Constraints
        
%         function [nonneg_infeas, nonneg_feas] = form_nonneg(obj, poly)
        function [nonneg_infeas, nonneg] = form_nonneg(obj, poly)
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
            nonneg.term = -vT; %vT <= 0
            
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
            
            nonneg_infeas = nonneg;            
%             nonneg_feas = nonneg;
%             nonneg_feas.init = v0 - poly.gamma;
%             nonneg_feas.term = t - vT; %minimum time control
        end    
        
        
        
        function cc_curr = make_psatz(obj, d, X, f, vars)
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
            
            cc_curr = coef_con(coeff, cons);
        end
        
        function cc_all = make_cons(obj, d, nonneg) %, feas)
            %make constraints for set disconnectedness program
            
%             if nargin < 4
%                 feas = 0;
%             end
            
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

            %occupation measure (+box, complement)
            if obj.options.scale
                Tsupp = t*(1-t);
            else
                Tsupp = t*(obj.options.Tmax - t);
            end          
           
            
            %terminal set            
            cc_term = coef_con();
%             if feas
%                 %minimum time control (minimize <t, mu_T>)
%                 %hopefully this is better conditioned than the feasibility
%                 %program
%                 X1_curr = X1_cell{i};
%                 X1_curr.ineq = [X1_curr.ineq; Tsupp];
%                 for i = 1:length(X1_cell)
%                     cc_term_curr = obj.make_psatz(d, X1_curr, nonneg.term, [t; x]);
%                     cc_term = [cc_term; cc_term_curr];
%                 end
%             else
                for i = 1:length(X1_cell)
                    cc_term_curr = obj.make_psatz(d, X1_cell{i}, nonneg.term, x);
                    cc_term = [cc_term; cc_term_curr];
                end
%             end
             
            
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
               
        
        
    end
end

