classdef lobe_plotter
    %LOBE_PLOTTER Visualization of Double-Lobe system
    %for set (dis)-connectedness
    
    properties
        FS_axis = 14;       %font size of axis
        FS_title = 16;      %font size of title
        
        f; %single set X
        
        out = [];
        out_sim = [];
        opt = [];
        func = struct('f', [], 'v', [], 'zeta', [], 'v0', [], 'v1', []);
        t;
        x;
        
        axlim = struct('x', [], 'y', [], 't', []);
        color0 = [0.4940, 0.1840, 0.5560];
        color1 = [0.4660, 0.6740, 0.1880];
    end
    
    methods
        function obj = lobe_plotter(opt, out, out_sim)
            %LOBE_PLOTTER Construct an instance of this class
            %   Detailed explanation goes here
            
            
            %properties of problem
            obj.opt = opt;
            
            %properties of solution
            if nargin >= 2
                obj.out = out;
            end
            
            if nargin >= 3
                obj.out_sim = out_sim;
            end
            
            %variables
            syms tv [1 1];
            syms xv [2 1];
            obj.t = tv;
            obj.x = xv;
            
            %symbolic evaluation only for the double-lobe
            %for plotting contours
            f_func = polyval_func(opt.X.ineq, opt.x);
            obj.func.X = f_func(xv);
            obj.func.v0= out.func.v0(xv);
            obj.func.v1= out.func.v1(xv);
            obj.func.v= out.func.v([tv; xv]);
            obj.func.zeta= out.func.zeta([tv; xv]);
            
            %plot axes
            obj.axlim.x= [-3, 2.5];
            obj.axlim.y= [-2, 2];
            if obj.opt.scale
                obj.axlim.t = [0, 1];
            else
                obj.axlim.t = [0, obj.opt.Tmax];
            end
        end
               
        
        function F = contour_2d(obj)
            %CONTOUR_2D Plot the double lobe and certificate in 2d
            %   Detailed explanation goes here
            F = figure(10);
            clf
            hold on
            
            %axis limits
            limits = [obj.axlim.x, obj.axlim.y];
            xlim(obj.axlim.x);
            ylim(obj.axlim.y);
    
            %sampled trajectories
            for i = 1:length(obj.out_sim)
            %     out_sim = set_walk(x0(), X_func, @() u_func(2), Tmax, dt);
                plot(obj.out_sim{i}.x(1, :), obj.out_sim{i}.x(2, :), 'c', 'HandleVisibility', 'off');
            end 
            
            %initial and final locations
            scatter(obj.opt.X0(1, :), obj.opt.X0(2, :), 100, 'ok', 'DisplayName', 'X0')
            scatter(obj.opt.X1(1, :), obj.opt.X1(2, :), 100, '*k', 'DisplayName', 'X1')

            fimplicit(obj.func.X == 0, limits,  'k', 'DisplayName','X')
            
            %contours of separation
            fimplicit(obj.func.v0 == 1, limits,'DisplayName','v(0, x) = 1', 'Color', obj.color0)
            fimplicit(obj.func.v1 == 0, limits,'DisplayName','v(T, x) = 0', 'Color', obj.color1)
    
            
            legend('location', 'northwest')
            
            title('Certificate of Disconnectedness', 'FontSize', obj.FS_title)
            xlabel('x_1')
            ylabel('x_2')
%             view(3)
            hold off                        
        end
        
        function F = contour_3d(obj)
            F = figure(11);
            clf
            hold on
            %axis limits
            
            limits = [obj.axlim.t, obj.axlim.x, obj.axlim.y];
            xlim(obj.axlim.t);
            ylim(obj.axlim.x);
            zlim(obj.axlim.y);

            
            %plot initial and final points
            n0 = size(obj.opt.X0, 2);
            n1 = size(obj.opt.X1, 2);
            scatter3(zeros(n0, 1), obj.opt.X0(1, :), obj.opt.X0(2, :), 100, 'ok', 'DisplayName', 'X0')
            scatter3(obj.axlim.t(2)*ones(n1, 1), obj.opt.X1(1, :), obj.opt.X1(2, :), 100, '*k', 'DisplayName', 'X1')

            for i = 1:n0
                plot3(obj.axlim.t, obj.opt.X0(1, i)*[1,1], obj.opt.X0(2, i)*[1,1], 'k', 'HandleVisibility', 'Off')
            end
            
            for i = 1:n1
                plot3(obj.axlim.t, obj.opt.X1(1, i)*[1,1], obj.opt.X1(2, i)*[1,1], 'k', 'HandleVisibility', 'Off')
            end
            
            fimplicit3(obj.func.v == 0, limits, 'MeshDensity', 120)
            
            
            xlabel('t')
            ylabel('x_1')
            zlabel('x_2')
            title('Certificate of Disconnectedness', 'FontSize', obj.FS_title)
        end
        
        %% separate these functions into an interface tomorrow
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
            
            Nzeta = length(obj.out.poly.zeta);
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
            Nzeta = length(obj.out.poly.zeta);
           
          
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
end

