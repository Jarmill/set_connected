classdef set_manager_partition < handle
    %SET_MANAGER_PARTITION 
    %Given sets X0 and X1 within a base space X,
    %determine if there exists a path between points in X0 and X1 (feas)
    %or there does not exist such a path (infeas)
    properties
        options;
        spacing; %number of locations at each dimension
        loc;     %the locations on the dimension
        limits;  %box limits for each cell
    end
    
    methods
        function obj = set_manager_partition(opt, spacing)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;

            obj.options = opt;
            if nargin < 2
                obj.spacing = ones(length(opt.x)+1,1);
            else
                if length(obj.spacing) == 1
                    obj.spacing = [obj.spacing; ones(length(opt.x), 1)]';
                else
                    obj.spacing = spacing;
                end
            end
            
%             obj.loc= cell(obj.spacing);
            
            if isempty(obj.options.t)
                obj.options.t = sdpvar(1, 1);
            end
            
            [obj.loc, obj.limits] = obj.make_locations(opt, spacing);
        end
        
        %% main execution routine
        
       
        function out = climb_connected(obj, order_range)
            %meta-algorithm: increase degree between d_range(1):d_range(2)
            
            
            if nargin == 1
                order_range = [1; 5];
            end
            
            if length(order_range) == 1
                order_range = [1; order_range];
            end
            
            for order_curr = order_range(1):order_range(2)
                out = obj.check_connected(2*order_curr);
                fprintf('order=%d, \t status=%s\n', order_curr, char(out.status))
                if (out.status ~= conn_status.Indeterminate)
                    break
                end
            end
            
        end
        
        function out = check_connected(obj, d)            
            prog_infeas = obj.make_program(d);
            out_infeas = solve_program(obj, prog_infeas);
            farkas= ~out_infeas.problem;
            if farkas
                status = conn_status.Disconnected;
            else
                status = conn_status.Indeterminate;          
            end
            
            out = out_infeas;
            out.status = status;
        end
        
        
        %% create the location cells of the partition
        
        function [loc, limits] = make_locations(obj, options, spacing)
            %create a grid of locations
            n = length(options.x);
            options.box = box_process(n, options.box);
            limits = obj.get_limits(options, spacing);
            spacing = reshape(spacing, 1, []);
            loc = cell(spacing);
            
            Ncell = prod(spacing);
            sub = cell(size(spacing));
            
            %create a location in every cell
            for i = 1:Ncell
                [sub{:}] = ind2sub(spacing, i);
%                 sublist = cell2num(sub);
                curr_time= limits{1}(sub{1} + [0, 1]);
                
%                 curr_box = cellfun(@(l, i) l{i}(sublist(i) + [0;1]), limits, num2cell([2:length(spacing)]));

                curr_box = zeros(size(options.box));
                for j = 1:n
                    curr_box(j, :) = limits{j+1}(sub{j+1} + [0, 1]);
                end                                                                
                
                
                loc{i} = set_location(options, curr_time, curr_box, i);
                
                %find the 'next' cell in each direction
                for j = 1:n+1
                    if sub{j} < spacing(j)
                        sub_new = sub;
                        sub_new{j} = sub{j} + 1;
                        ind_next = sub2ind(spacing, sub_new{:});
                        
                        loc{i}.next{j} = ind_next;
                    end
                end
            end                
        end
        
        function limits = get_limits(obj, options, spacing)
            %form the limits (ranges) of each location
%             loc = cell(spacing);
            
            limits = cell(length(spacing), 1);
            limits{1} = linspace(0, options.Tmax, spacing(1)+1);
            for i = 2:length(spacing)
                limits{i} = linspace(options.box(i-1, 1), options.box(i-1, 2), spacing(i)+1);
            end
        end
            
        %% generate adjacency constraints between cells
        function [con_adj, coeff_adj] = make_adjacency_con(obj, d, poly_var,loc)
            % inputs: poly_var and loc are cells. The attribute loc.poly
            % may not yet have been bound to poly_var
            
            con_adj = [];
            
            coeff_adj = [];
            
            t = obj.options.t;
            x = obj.options.x;
            vars = [t; x];
            
            spacing = size(loc);
            Ncell = prod(spacing);
            
            for i = 1:Ncell
                loc_curr = loc{i};
                for j = 1:length(spacing)
                    ind_next = loc_curr.next{j};
                    %if there is an adjacent face, form the equality
                    %constraints between polynomials
                    if ~isempty(ind_next)
                       
                        %polynomials at current cell
                        poly_curr = loc_curr.get_adjacency_poly(poly_var{i}, j-1, 1);
                        
                        %polynomials at next cell
                        loc_next = loc{ind_next};                        
                        poly_next = loc_next.get_adjacency_poly(poly_var{ind_next}, j-1, 0);
                        
                        if j == 1
                            %time cell
                            X_shared = loc_curr.get_X();
                            
                            %v should increase along the time transition
                            pos_time_jump = poly_next - poly_curr;
                            cons_time = [];
                            coeff_time = [];
%                             [con_u_curr, coeff_u_curr] = obj.make_psatz(d, X_curr, X_shared, [x]);
                            for i = 1:length(X_shared)
                                [p_time_curr, con_time_curr, coeff_time_curr] = constraint_psatz(pos_time_jump, X_shared{i}, x, d); 
                                cons_time = [cons_time; con_time_curr];
                                coeff_time = [coeff_time; coeff_time_curr];
                            end
%                     con_u_curr = con_u_curr:['Cell ', num2str(obj.id), ' Lie input ', num2str(j)];
%                     con_u = [con_u; con_u_curr];
%                     coeff_u = [coeff_u; coeff_u_curr];
                            con_adj_curr = cons_time;
                            coeff_adj = [coeff_adj; coeff_time_curr];

                        else
                            %space cell
                            coeff_pv = coefficients(poly_curr - poly_next, vars);
                            con_adj_curr = coeff_pv == 0;
                        end
                        con_adj = [con_adj; con_adj_curr:['Cell ', num2str(i), '-', num2str(ind_next)]];
%                         end
                    end
                end
                
            end
        end        
        
        
%         function [con_adj] = adjacency_space(obj, loc_curr, loc_next)
%             %adjacency_space constraints ensuring that v is the same
%             %between space cells. Unfortunately equalities are required,
%             %there are no inequality relaxations here
%             
%             v_curr = loc_vu
%             
%         end
        
        %% form and solve the program
        function [prog_mgr] = make_program(obj, d)
            %combine together constraints from all locations 
            
%             prog_mgr = [];
            
            Ncell = prod(size(obj.loc));
%             poly_var = 
            cons_mgr = [];
            coeff_mgr = [];
            
%             for i = 1:Ncell
%                 prog_loc = obj.loc{i}r          
            
            prog_all = cellfun(@(l) l.make_program(d), obj.loc, 'uniformoutput', false);
            poly_var = cellfun(@(p) p.poly, prog_all, 'uniformoutput', false);
            nonneg_var = cellfun(@(p) p.nonneg, prog_all, 'uniformoutput', false);
            
            con_all = cellfun(@(p) p.cons, prog_all, 'uniformoutput', false);
            con_all = [con_all{:}];
            
            coeff_all = cellfun(@(p) p.coeff, prog_all, 'uniformoutput', false);
            coeff_all = reshape([coeff_all{:}], [], 1);
            
            
            [con_adj, coeff_adj] = obj.make_adjacency_con(d, poly_var,obj.loc);
            prog_mgr = struct('objective', 0, 'coeff', [coeff_all; coeff_adj], 'cons', [con_all; con_adj]);
%             prog_mgr = struct('nonneg', nonneg_var, 'poly', poly_var, 'cons', ...
%                 [con_all; con_adj], 'coeff', coeff_all, 'objective', 0);

            prog_mgr.poly = poly_var; 
            prog_mgr.nonneg = nonneg_var;
        end
        
        function [out] = solve_program(obj, prog)
            %solve SOS program in YALMIP, return solution    
            
            opts = sdpsettings('solver', obj.options.solver, 'verbose', obj.options.verbose);
            opts.sos.model = 2;
            
            [sol, monom, Gram, residual] = solvesos(prog.cons, prog.objective, opts, prog.coeff);
            
            out = struct('poly', [], 'problem', sol.problem, 'sol', [], 'block', [], 'func', [], ...
                'n', length(obj.options.x), 'limits', []); 
            
            %putting a cell in a struct constructor will yield a struct
            %array. This is undesirable behavior
            out.limits = obj.limits;
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
    
%         function
        
        function [poly_rec, func_rec] = recover_poly(obj, poly, nonneg)
            %distribute out recovery to the locations
            
            Ncell = prod(size(obj.loc));
%             poly_rec = cell(size(obj.loc));
%             func_rec = cell(size(obj.loc));

            rec_ind = @(i) obj.loc{i}.recover_poly(poly{i}, nonneg{i});
            [poly_rec, func_cell] = arrayfun(rec_ind, ...
                1:Ncell, 'UniformOutput', false);
            
%             poly_rec = cellfun(@(o) o{1}, out_rec{:});                       
%             func_cell = cellfun(@(o) o{2}, out_rec{:});
            
            %export the function evaluators
            func_rec = struct;
            
            func_rec.func_cell = func_cell;
            
            %determine the grid spacing
            spacing = cellfun(@length, obj.limits)-1;
            dlim = cellfun(@(l) l(2) - l(1), obj.limits);
            start = cellfun(@(l) l(1), obj.limits);
            
            func_rec.grid_handle = @(x) grid_ind(x, dlim, start, spacing);
      
            
            if obj.options.scale
                scale_weight = obj.options.Tmax;
            else
                scale_weight = 1;
            end
            
            %handles to evaluate the polynomials
            
            func_rec.v      = @(vars) func_rec.func_cell{func_rec.grid_handle([vars(1)/scale_weight; vars(2:end)])}.v([vars(1)/scale_weight; vars(2:end)]);
            func_rec.v0     = @(vars) func_rec.func_cell{func_rec.grid_handle([zeros(1,size(vars,2));vars])}.v0(vars);
            func_rec.v1     = @(vars) func_rec.func_cell{func_rec.grid_handle([obj.options.Tmax*ones(1,size(vars,2));vars])}.v1(vars);
            func_rec.zeta   = @(vars) func_rec.func_cell{func_rec.grid_handle([vars(1)/scale_weight; vars(2:end)])}.zeta([vars(1)/scale_weight; vars(2:end)]);
            func_rec.nonneg = @(vars) func_rec.func_cell{func_rec.grid_handle([vars(1)/scale_weight; vars(2:end)])}.nonneg([vars(1)/scale_weight; vars(2:end)]);
            
            %convert to a cellfun
%             for i = 1:Ncell
%                 [poly_rec{i}, func_rec{i}] = obj.loc{i}.recover_poly(poly{i}, nonneg{i});
%             end
            
        end
    end
end

