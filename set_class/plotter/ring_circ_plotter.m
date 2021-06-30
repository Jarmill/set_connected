classdef ring_circ_plotter < set_plotter_interface
    %RING_CIRC plotter visualization of ring-circle plotter for double-lobe
    %system
    
    properties        
        axlim = struct('x', [-1.5,1.5], 'y', [-1.5,1.5], 't', []);

        rad0 = [];
        circle; %description of set X
    end
    
    methods
        function obj = ring_circ_plotter(opt, out, out_sim)
            %LOBE_PLOTTER Construct an instance of this class
            %   Detailed explanation goes here
            
             
            %properties of solution
            if nargin < 2
                out = [];
            end
            
            if nargin < 3
                out_sim = [];
            end
            
            obj = obj@set_plotter_interface(opt,out, out_sim);
            %properties of problem

            
            
            %variables
%             syms tv [1 1];
  
            %plot axes
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
%             fimplicit(obj.func.X == 0, limits,  'k', 'DisplayName','X')
            
            th_fine = linspace(0, 2*pi, 200);
            circ_fine = [cos(th_fine); sin(th_fine)];
    

            color_int = 0.8*[1,1,1];
            patch(obj.circle.R_ring_outer*circ_fine(1, :), obj.circle.R_ring_outer*circ_fine(2, :), color_int, 'LineWidth', 3, 'DisplayName', 'X')
            patch(obj.circle.R_ring_inner*circ_fine(1, :), obj.circle.R_ring_inner*circ_fine(2, :), 'w', 'LineWidth', 3, 'HandleVisibility', 'off')
            patch(obj.circle.R_circ*circ_fine(1, :), obj.circle.R_circ*circ_fine(2, :), color_int, 'LineWidth', 3, 'HandleVisibility', 'off')

            if isnumeric(obj.opt.X0)
                if size(obj.opt.X0, 1)==1
                    plot(obj.opt.X0*circ_fine(1, :), obj.opt.X0*circ_fine(2, :), 'color', obj.color0, 'LineWidth', 3, 'DisplayName', 'X0')
                else
                    scatter(obj.opt.X0(1, :), obj.opt.X0(2, :), 100, 'ok', 'DisplayName', 'X0')
                end
            end
            
            if isnumeric(obj.opt.X1)
                if size(obj.opt.X0, 1)==1
                    plot(obj.opt.X1*circ_fine(1, :), obj.opt.X1*circ_fine(2, :), 'color', obj.color1, 'LineWidth', 3, 'DisplayName', 'X1')
                else
                    scatter(obj.opt.X1(1, :), obj.opt.X1(2, :), 100, '*k', 'DisplayName', 'X1')
                end
            end
            
            
            title('Ring Circle set to analyze', 'FontSize', obj.FS_title)
            
            xlabel('x_1')
            ylabel('x_2')
            axis square
            legend('location', 'northwest')
            
        end
        
        function F = contour_2d(obj)
            %CONTOUR_2D Plot the double lobe and certificate in 2d
            %   Detailed explanation goes here
            [F, a] = obj.set_plot();
            hold on
            
            %axis limits
            limits = [obj.axlim.x, obj.axlim.y];
    
            %sampled trajectories
            for i = 1:length(obj.out_sim)
            %     out_sim = set_walk(x0(), X_func, @() u_func(2), Tmax, dt);
                plot(a,obj.out_sim{i}.x(1, :), obj.out_sim{i}.x(2, :), 'c', 'HandleVisibility', 'off');
            end 
            
            fcell = obj.out.func.func_cell;
            for i = 1:length(fcell)
                curr_v = fcell{i}.v;
                curr_lim = reshape(fcell{i}.box', 1, []);
%                 fsurf(@(t,r) curr_v([t/opt.Tmax;r]), curr_lim,'DisplayName','v(t, x)')
                fcontour(@(t,r) curr_v([t/obj.opt.Tmax;r]), curr_lim, 'k', 'LineWidth', 4, 'DisplayName','v(t, x)=0', 'LevelList', 0)
            end
            
            %initial and final locations
%             scatter(a,obj.opt.X0(1, :), obj.opt.X0(2, :), 100, 'ok', 'DisplayName', 'X0')
%             scatter(a,obj.opt.X1(1, :), obj.opt.X1(2, :), 100, '*k', 'DisplayName', 'X1')

            
            %contours of separation
%             fimplicit(a,obj.func.v0 == opt.epsilon, limits,'DisplayName','v(0, x) = 1', 'Color', obj.color0)
%             fimplicit(a,obj.func.v1 == 0, limits,'DisplayName','v(T, x) = 0', 'Color', obj.color1)
    
            
            
            
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
            
%             fimplicit3(obj.func.v == 0, limits, 'MeshDensity', 120)
            
            
            xlabel('t')
            ylabel('x_1')
            zlabel('x_2')
            title('Certificate of Disconnectedness', 'FontSize', obj.FS_title)
        end
       
    end
end

