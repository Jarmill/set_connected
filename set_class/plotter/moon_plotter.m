classdef moon_plotter < set_plotter_interface
    %MOON_PLOTTER Visualization of bicuspid system with circles
    %for set (dis)-connectedness
    
    properties
        
        %moon contour
        inner_x = []; %(x - inner_x)^2 + y^2 >= inner_rad^2 (assume inner_x is positive)
        inner_rad = [];
        
        %other circles
        circle = struct('center', [], 'rad', []);
    end
    
    methods
        function obj = moon_plotter(opt, out, out_sim)
            %MOON_PLOTTER Construct an instance of this class
            %   Detailed explanation goes here
            
            
            if nargin < 2
                out = [];
            end
            
            if nargin < 3
                out_sim = [];
            end
            
            obj = obj@set_plotter_interface(opt,out, out_sim);
            %properties of problem          
           
            
            
            %plot axes
            obj.axlim.x= [-3, 2.5];
            obj.axlim.y= [-2, 2];
            if obj.opt.scale
                obj.axlim.t = [0, 1];
            else
                obj.axlim.t = [0, obj.opt.Tmax];
            end
            
            if out.problem
                return
            end
        end
               
        
        function [x_circ, x_moon] = moon_profile(obj, inner_center, inner_rad, N)
            %draw points on the moon
            if nargin < 4
                N = 300;
            end
            theta = linspace(0, 2*pi, N);
            x_circ = [cos(theta); sin(theta)];
            
            x_outer = x_circ;
            x_inner = inner_rad*x_circ + [inner_center; 0];
            
%             f_outer = 1- sum(x_circ(1, :).^2 + x_circ(2, :).^2);
            f_outer = -inner_rad^2+ (x_outer(1, :)-inner_center).^2 + x_outer(2, :).^2;
            
            f_inner = 1- x_inner(1, :).^2 - x_inner(2, :).^2;
            
            x_outer_filter = x_outer(:, f_outer >= 0);
            x_inner_filter = x_inner(:, f_inner >= 0);
            
            x_moon = [x_outer_filter, x_inner_filter(:, end:-1:1), x_outer_filter(:, 1)];                                    
        end
        
        function [F, a] = set_plot(obj)
            F = figure(9);
            clf
            a = axes;
            hold on
            
            
            limits = [obj.axlim.x, obj.axlim.y];

            [x_circ, x_moon] = obj.moon_profile(obj.inner_x, obj.inner_rad);
            plot(x_moon(1, :), x_moon(2, :), 'k', 'DisplayName', 'X')
            
            for i = 1:length(obj.circle.rad)
                rx = obj.circle.rad(i);
                cx = obj.circle.center(:, i);
                x_circ_curr = rx*x_circ + cx;
                plot(x_circ_curr(1, :), x_circ_curr(2, :), 'k', 'HandleVisibility', 'off')
            end
%             fimplicit(obj.func.X == 0, limits,  'k', 'DisplayName','X')
            
            %circles go here
            
            scatter(obj.opt.X0(1, :), obj.opt.X0(2, :), 100, 'ok', 'DisplayName', 'X0')
            scatter(obj.opt.X1(1, :), obj.opt.X1(2, :), 100, '*k', 'DisplayName', 'X1')

            
            title('Moon/Circle set to analyze', 'FontSize', obj.FS_title)
            xlabel('x_1')
            ylabel('x_2')
            xlim(obj.axlim.x);
            ylim(obj.axlim.y);
%             view(3)
            hold off        
        end
        
        function F = contour_2d(obj)
            %CONTOUR_2D Plot the double lobe and certificate in 2d
            %   Detailed explanation goes here
%             F = figure(10);
%             clf
%             hold on

            [F, a] = obj.set_plot();
%             a = axes;
            %axis limits
            
            limits = [obj.axlim.x, obj.axlim.y];
            
            %contours of separation
            hold on
            fimplicit(a, obj.func.v0 == 1, limits,'DisplayName','v(0, x) = 1', 'Color', obj.color0)
            fimplicit(a, obj.func.v1 == 0, limits,'DisplayName','v(T, x) = 0', 'Color', obj.color1)
    
            
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
            
%             MD = 120;
            fimplicit3(obj.func.v == 0, limits, 'MeshDensity', 120)
            
            
            xlabel('t')
            ylabel('x_1')
            zlabel('x_2')
            title('Certificate of Disconnectedness', 'FontSize', obj.FS_title)
        end
    end
end

