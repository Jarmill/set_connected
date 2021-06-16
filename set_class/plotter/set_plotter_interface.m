classdef set_plotter_interface < handle
    %set_plotter_interface provides functions to visualize certificates of
    %set connectedness
    
    properties
        FS_axis = 14;       %font size of axis
        FS_title = 16;      %font size of title1
        
        
    end
    
    properties(Access=protected)
        opt;            %options for plotting
        out;            %result of optimization
        out_sim;        %trajectories
        
        vars=struct('t',[],'x',[]);           %symbolic variables for plotting
%         func = struct('v', [], 'zeta', [], 'v0', [], 'v1', []);   %functions in symbolic variables
        
        X_func;

        %colors for v contours at times 0 and T
        color0 = [0.4940, 0.1840, 0.5560];
        color1 = [0.4660, 0.6740, 0.1880];
    end
    
    methods
        function obj = set_plotter_interface(opt,out, out_sim)
            %fill in data for the certificate and sampled trajectories
            obj.opt = opt;
            obj.out = out;
            obj.out_sim = out_sim;
            
            
            %set up symbolic variables
%             n = out.n;
            
           
%             syms tv [1 1];
%             syms xv [n 1];
%             obj.vars.t = tv;
%             obj.vars.x = xv;
            
%             if ~isempty(obj.out.func)
%                 obj.func.v0= out.func.v0(xv);
%                 obj.func.v1= out.func.v1(xv);
%                 obj.func.v= out.func.v([tv; xv]);
%                 obj.func.zeta= out.func.zeta([tv; xv]);
%             end      
            
            
        end
        
        function F = state_plot(obj)
            F = figure(19);
            clf
            nx = size(obj.out_sim{1}.x, 2);
            nplt = nx;
            ax = cell(nplt, 1);
            tiledlayout(nplt,1);
            for k = 1:nplt
                ax{k} = nexttile;
                title(['State ', num2str(k), ' vs. Time'], 'FontSize', obj.FS_title);   
                ax_loc_curr = ['$x_', num2str(k), '$'];
                xlabel('time', 'FontSize', obj.FS_axis)
                ylabel(ax_loc_curr , 'interpreter', 'latex', 'FontSize', obj.FS_axis);
                hold on
            end
            
            for i = 1:length(obj.out_sim)
                tcurr = obj.out_sim{i}.t;
                xcurr = obj.out_sim{i}.x;
                for k = 1:nplt
                   plot(ax{k}, tcurr, xcurr(:, k), 'c') 
                end
            end
            
        end
        
        function [F, limits] = state_plot_2(obj, box_lim)
            F = figure(30);
            clf
            hold on 
            for i = 1:length(obj.out_sim)
                tcurr = obj.out_sim{i}.t;
                xcurr = obj.out_sim{i}.x;
                plot(xcurr(:, 1), xcurr(:, 2), 'c')
            end
            
            %don't mark the initial conditions
%             for i = 1:length(obj.out_sim)
%                 tcurr = obj.out_sim{i}.t;
%                 xcurr = obj.out_sim{i}.x;
%                 scatter(xcurr(1, 1), xcurr(1, 2), 100, 'k')
%             end                        
            
            xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Phase Plane', 'FontSize', obj.FS_title);   
            
            
            if (nargin == 2 ) && ~isempty(box_lim)
                box = box_process(2, box_lim);
                xlim(box(1, :));
                ylim(box(2, :));
                limits = [box(1, :), box(2, :)];
            else
                limits = [xlim, ylim];
            end
            
            pbaspect([diff(xlim), diff(ylim), 1])
            
        end
 
        function [F, limits] = state_plot_3(obj, box_lim)
            F = figure(30);
            clf
            hold on 
            for i = 1:length(obj.out_sim)
                tcurr = obj.out_sim{i}.t;
                xcurr = obj.out_sim{i}.x;
                plot3(xcurr(:, 1), xcurr(:, 2), xcurr(:, 3), 'c')
            end
            
            %don't mark the initial conditions
%             for i = 1:length(obj.out_sim)
%                 tcurr = obj.out_sim{i}.t;
%                 xcurr = obj.out_sim{i}.x;
%                 scatter(xcurr(1, 1), xcurr(1, 2), 100, 'k')
%             end                        
            
            xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Phase Plane', 'FontSize', obj.FS_title);   
            
            
            if (nargin == 2 ) && ~isempty(box_lim)
                box = box_process(3, box_lim);
                xlim(box(1, :));
                ylim(box(2, :));
                zlim(box(3, :));
                limits = [box(1, :), box(2, :), box(3, :)];
            else
                limits = [xlim, ylim, zlim];
            end
            
            pbaspect([diff(xlim), diff(ylim), diff(zlim)])
            view(3);
        end
        
        function F = v_plot(obj)

            F = figure(18);
            clf
            subplot(2, 1, 1)
            hold on

            for j = 1:length(obj.out_sim)
                osc = obj.out_sim{j};                    
                plot(osc.t_traj, osc.v, 'c');

            end
            
%             plot(xlim, [1;1]*obj.out.poly.gamma, '--r', 'LineWidth', 3)
            xlabel('time', 'FontSize', obj.FS_axis)
            ylabel('$v(t,x)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Auxiliary Function', 'FontSize', obj.FS_title);   
            
            subplot(2,1,2)
            hold on
            for j = 1:length(obj.out_sim)
                osc = obj.out_sim{j};                    
                plot(osc.t_traj(1:end-1), diff(osc.v), 'c');

            end
            xlabel('time', 'FontSize', obj.FS_axis)
            ylabel('$\Delta v(t,x)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Auxiliary Function Increase', 'FontSize', obj.FS_title);   
            plot(xlim, [0,0], 'k:') 
            
                end
        
        function F = nonneg_zeta(obj)
            %plot the nonnegative slack functions zeta
            
            Nzeta = obj.out.n;
            if Nzeta
            F = figure(21);
            clf
            
            
            hold on
            for j = 1:length(obj.out_sim)
                osc = obj.out_sim{j};                    
                plot(osc.t_traj, osc.nonneg(end-Nzeta:end, :), 'c')
            end

            ax_loc_curr = ['$\zeta(t,x)$'];
            xlabel('time', 'FontSize', obj.FS_axis)
            ylabel(ax_loc_curr , 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title(['Constraint Slack'], 'FontSize', obj.FS_title);   
            plot(xlim, [0,0], 'k:')
            else
                F = [];
            end
            
        end
        
        function F = nonneg_traj(obj)

            F = figure(20);
            clf
%             N0 = length(obj.manager.opts.get_X_init());
            
            %no other switching present, so there is only one system for
            %the occupation measure
            Nzeta = obj.out.n;
           
          
            ax_title = {'Nonpositive v (initial)', 'Reachability Indicator (terminal)', 'Reachability Indicator (terminal slack)', 'Decreasing v (occupation)'};
            
            tiledlayout(Nzeta+1, 1);
            for i = 1:Nzeta+1
                nexttile
                hold on
                for j = 1:length(obj.out_sim)
                    osc = obj.out_sim{j};                    
                    plot(osc.t_traj, osc.nonneg(i, :), 'c');
                end
                
                xlabel('time', 'FontSize', obj.FS_axis)
                plot(xlim, [0,0], 'k:')
                
                if i==1
                    ylabel(['$\partial_t v - \nabla_x \cdot v -  1 \cdot \zeta_i$'], 'interpreter', 'latex', 'FontSize', obj.FS_axis);
                    title(['Decreasing v'],'FontSize', obj.FS_title);   
                else
                    si = num2str(i-1);
                    ylabel(['$\zeta_', si, '+2 \nabla_x \cdot v$'], 'interpreter', 'latex', 'FontSize', obj.FS_axis);
                    title(['Input ', si, ' Bounding'],'FontSize', obj.FS_title);   
                end
            end
         
        end
    
    end
    
    methods(Abstract) 
        %draw the set X;
        set_plot(obj, F);
    end
end

