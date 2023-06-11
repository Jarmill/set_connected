classdef elliptic_curve_plotter
    %ELLIPTIC_CURVE_PLOTTER Visualization of bicuspid system with circles
    %for set (dis)-connectedness
    
    properties
        FS_axis = 14;       %font size of axis
        FS_title = 16;      %font size of title
        
        %Weierstrass-form elliptic curve
        %y^2 = x^3 + a x^2 + b
        a =0.8; %parameter of bicuspid
        b =-0.05; %parameter of bicuspid
        %other circles
        circle_center = [];
        circle_rad = [];
        
        out = [];
        opt = [];
        func = struct('f', [], 'v', [], 'zeta', [], 'v0', [], 'v1', []);
        t;
        x;
        
        axlim = struct('x', [], 'y', [], 't', []);
        color0 = [0.4940, 0.1840, 0.5560];
        color1 = [0.4660, 0.6740, 0.1880];
    end
    
    methods
        function obj = elliptic_curve_plotter(opt, out)
            %ELLIPTIC_CURVE_PLOTTER Construct an instance of this class
            %   Detailed explanation goes here
            
            
            %properties of problem
            obj.opt = opt;
            
            %properties of solution
            if nargin == 2
                obj.out = out;
            end
            
            %variables
            syms tv [1 1];
            syms xv [2 1];
            obj.t = tv;
            obj.x = xv;
            
            %symbolic evaluation only for the double-lobe
            %for plotting contours
%             f_func = polyval_func(opt.X.ineq, opt.x);
%             obj.func.X = f_func(xv);
            %bicuspid
            f_func = @(x, y) y.^2 - x.^3 - obj.a*x.^2 - obj.b;
            obj.func.X = f_func(xv(1), xv(2));
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
               
        
        function F = set_plot(obj)
            F = figure(9);
            clf
            hold on
            
            limits = [obj.axlim.x, obj.axlim.y];

     
            
            fimplicit(obj.func.X == 0, limits,  'k', 'DisplayName','X')
            
            %circles go here
            
            title('Elliptic Curve set to analyze', 'FontSize', obj.FS_title)
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
            F = figure(10);
            clf
            hold on
            
            %axis limits
            limits = [obj.axlim.x, obj.axlim.y];
            xlim(obj.axlim.x);
            ylim(obj.axlim.y);
    
            %initial and final locations
            %only if points (correct this later)
            scatter(obj.opt.X0(1, :), obj.opt.X0(2, :), 100, 'ok', 'DisplayName', 'X0')
            scatter(obj.opt.X1(1, :), obj.opt.X1(2, :), 100, '*k', 'DisplayName', 'X1')

            fimplicit(obj.func.X == 0, limits,  'k', 'DisplayName','X', 'edgecolor', 'none')
            
            %contours of separation
            fimplicit(obj.func.v0 == 1, limits,'DisplayName','v(0, x) = 1', 'Color', obj.color0)
%             fimplicit(obj.func.v0 == 0, '--', limits,'DisplayName','v(0, x) = 0', 'Color', obj.color0)
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
            scatter3(zeros(n0, 1), obj.opt.X0(1, :), obj.opt.X0(2, :), 200, 'or', 'DisplayName', 'X0')
            scatter3(obj.axlim.t(2)*ones(n1, 1), obj.opt.X1(1, :), obj.opt.X1(2, :), 200, '*r', 'DisplayName', 'X1')

            for i = 1:n0
                plot3(obj.axlim.t, obj.opt.X0(1, i)*[1,1], obj.opt.X0(2, i)*[1,1], 'k', 'HandleVisibility', 'Off')
            end
            
            for i = 1:n1
                plot3(obj.axlim.t, obj.opt.X1(1, i)*[1,1], obj.opt.X1(2, i)*[1,1], 'k', 'HandleVisibility', 'Off')
            end
            
            fimplicit3(obj.func.v == 0, limits, 'MeshDensity', 160)
            
            
            xlabel('t')
            ylabel('x_1')
            zlabel('x_2')
            title('Certificate of Disconnectedness', 'FontSize', obj.FS_title)
        end
    end
end

