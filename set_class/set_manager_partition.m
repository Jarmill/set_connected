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
        function con_adj = make_adjacency_con(obj, poly_var,loc)
            % inputs: poly_var and loc are cells. The attribute loc.poly
            % may not yet have been bound to poly_var
            
            con_adj = [];
            
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
                        
                        coeff_pv = coefficients(poly_curr - poly_next, vars);
                        con_adj_curr = coeff_pv == 0;
                        
                        con_adj = [con_adj; con_adj_curr:['Cell ', num2str(i), '-', num2str(ind_next)]];
                    end
                end
                
            end
        end        
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
            coeff_all = [coeff_all{:}];
            
            
            con_adj = obj.make_adjacency_con(poly_var,obj.loc);
            prog_mgr = struct('objective', 0, 'coeff', coeff_all, 'cons', [con_all; con_adj]);
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
            
            out = struct('poly', [], 'problem', sol.problem, 'sol', [], 'block', [], 'func', [], 'n', length(obj.options.x));
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
    
        function [poly_rec, func_rec] = recover_poly(obj, poly, nonneg)
            %distribute out recovery to the locations
            
            Ncell = prod(size(obj.loc));
            poly_rec = cell(size(obj.loc));
            func_rec = cell(size(obj.loc));
            
            %convert to a cellfun
            for i = 1:Ncell
                [poly_rec{i}, func_rec{i}] = obj.loc{i}.recover_poly(poly{i}, nonneg{i});
            end
            
        end
    end
end

