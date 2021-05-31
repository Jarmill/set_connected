classdef lobe_plotter < set_plotter_interface
    %LOBE_PLOTTER Visualization of Double-Lobe system
    %for set (dis)-connectedness
    
    properties
        opt = [];
        
        axlim = struct('x', [], 'y', [], 't', []);
        color0 = [0.4940, 0.1840, 0.5560];
        color1 = [0.4660, 0.6740, 0.1880];
    end
    
    methods
        function obj = lobe_plotter(opt, out, out_sim)
            %LOBE_PLOTTER Construct an instance of this class
            %   Detailed explanation goes here
            
             
            %properties of solution
            if nargin < 2
                out = [];
            end
            
            if nargin < 3
                out_sim = [];
            end
            
            obj = obj@set_plotter_interface(opt, out, out_sim);
            %properties of problem
%             obj.opt = opt;
            

            
            %symbolic evaluation only for the double-lobe
            %for plotting contours
            f_func = polyval_func(opt.X.ineq, opt.x);
            obj.func.X = f_func(obj.vars.x);
  
            %plot axes
            obj.axlim.x= [-3, 2.5];
            obj.axlim.y= [-2, 2];
            if obj.opt.scale
                obj.axlim.t = [0, 1];
            else
                obj.axlim.t = [0, obj.opt.Tmax];
            end
        end
               
         
        function [F, a] = set_plot(obj)
            F = figure(10);
            clf
            
            a = axes;
            hold on
            limits = [obj.axlim.x, obj.axlim.y];
            xlim(obj.axlim.x);
            ylim(obj.axlim.y);
            fimplicit(obj.func.X == 0, limits,  'k', 'DisplayName','X')
            
            title('Double Lobe set to analyze', 'FontSize', obj.FS_title)
            
            xlabel('x_1')
            ylabel('x_2')
            
        end
        
        function F = contour_2d(obj)
            %CONTOUR_2D Plot the double lobe and certificate in 2d
            %   Detailed explanation goes here
            [F, a] = obj.set_plot;
            hold on
            
            %axis limits
             limits = [obj.axlim.x, obj.axlim.y];
    
            %sampled trajectories
            for i = 1:length(obj.out_sim)
            %     out_sim = set_walk(x0(), X_func, @() u_func(2), Tmax, dt);
                plot(a,obj.out_sim{i}.x(1, :), obj.out_sim{i}.x(2, :), 'c', 'HandleVisibility', 'off');
            end 
            
            %initial and final locations
            scatter(a,obj.opt.X0(1, :), obj.opt.X0(2, :), 100, 'ok', 'DisplayName', 'X0')
            scatter(a,obj.opt.X1(1, :), obj.opt.X1(2, :), 100, '*k', 'DisplayName', 'X1')

            
            %contours of separation
            fimplicit(a,obj.func.v0 == 1, limits,'DisplayName','v(0, x) = 1', 'Color', obj.color0)
            fimplicit(a,obj.func.v1 == 0, limits,'DisplayName','v(T, x) = 0', 'Color', obj.color1)
    
            
            legend('location', 'northwest')
            
            title('Certificate of Disconnectedness', 'FontSize', obj.FS_title)

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
       
    end
end

