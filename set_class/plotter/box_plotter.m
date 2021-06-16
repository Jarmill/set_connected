classdef box_plotter < set_plotter_interface
    %BOX_PLOTTER Visualization of assemblage-of-boxes system
    %for set (dis)-connectedness
    
    properties
%         opt = [];

        boxes = {};
        
        axlim = struct('x', [], 'y', [], 't', []);

    end
    
    methods
        function obj = box_plotter(opt, out, out_sim)
            %box_PLOTTER Construct an instance of this class
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
%             f_func = polyval_func(opt.X.ineq, opt.x);

            obj.axlim.x = obj.opt.box(1, :);
            obj.axlim.y = obj.opt.box(2, :);
            
%             if obj.opt.scale
%                 obj.axlim.t = [0, 1];
%             else
                obj.axlim.t = [0, obj.opt.Tmax];
%             end
        end
               
         
        function [F, a] = set_plot(obj)
            F = figure(10);
            clf
            
            a = axes;
            hold on
            limits = [obj.axlim.x, obj.axlim.y];
            xlim(obj.axlim.x);
            ylim(obj.axlim.y);
%             fimplicit(@(x,y) obj.X_func([x;y]), limits,  'k', 'DisplayName','X')
 
            color_int = 0.8*[1,1,1];

            %plot the boxes
            for i = 1:length(obj.boxes)
                box_curr = obj.boxes{i};
                
                x_curr = box_curr(1,[1,2,2,1,1]);
                y_curr = box_curr(2,[1,1,2,2,1]);
                if i == 1 
                    patch(x_curr, y_curr, color_int, 'edgecolor', 'none', 'DisplayName', 'X');
                else
                    patch(x_curr, y_curr, color_int, 'edgecolor', 'none', 'HandleVisibility', 'off');
                end
            end
            
            title('Boxes to analyze', 'FontSize', obj.FS_title)
            
            xlabel('x_1')
            ylabel('x_2')
            
        end
        
        function [F,a] = traj_2d(obj)
            [F, a] = obj.set_plot;
            hold on
           
            %sampled trajectories
            for i = 1:length(obj.out_sim)
                plot(a,obj.out_sim{i}.x(1, :), obj.out_sim{i}.x(2, :), 'c', 'HandleVisibility', 'off');
            end 
            
            %initial and final locations
            scatter(a,obj.opt.X0(1, :), obj.opt.X0(2, :), 100, 'ok', 'DisplayName', 'X0')
            scatter(a,obj.opt.X1(1, :), obj.opt.X1(2, :), 100, '*k', 'DisplayName', 'X1')

        end
        
        function F = contour_2d(obj, tdiv)
            %CONTOUR_2D Plot the double lobe and certificate in 2d
            %   Detailed explanation goes here
            [F, a] = obj.set_plot;
                        %initial and final locations
            scatter(a,obj.opt.X0(1, :), obj.opt.X0(2, :), 100, 'ok', 'DisplayName', 'X0')
            scatter(a,obj.opt.X1(1, :), obj.opt.X1(2, :), 100, '*k', 'DisplayName', 'X1')

            hold on
            
            %axis limits
            limits = [obj.axlim.x, obj.axlim.y];
    
            
            %contours of separation
            
            v0_sep = @(x,y) obj.out.func.v0([x;y])-obj.opt.epsilon;
            v0_sep_vec = @(x,y)cell2mat(arrayfun(@(i)v0_sep(x(i), y(i)),...
    (1:length(x)),'UniformOutput',false));
            fimplicit(a, v0_sep_vec, limits,'DisplayName',['v(0, x) = ', num2str(obj.opt.epsilon)], 'Color', obj.color0)
            
            v1_sep = @(x,y) obj.out.func.v1([x;y]);
            v1_sep_vec = @(x,y)cell2mat(arrayfun(@(i)v1_sep(x(i), y(i)),...
    (1:length(x)),'UniformOutput',false));
            fimplicit(a,v1_sep_vec, limits,'DisplayName','v(T, x) = 0', 'Color', obj.color1)
    
            if nargin == 2
                %additional contours
                tspan = linspace(obj.axlim.t(1), obj.axlim.t(2), 2+tdiv);
                for i = 1:tdiv
                    tcurr = tspan(i+1);
                    v_sep = @(x,y) obj.out.func.v([tcurr;x;y]);
                    v_sep_vec = @(x,y)cell2mat(arrayfun(@(i)v_sep(x(i), y(i)),...
                    (1:length(x)),'UniformOutput',false));
                    vcurr = fimplicit(a,v_sep_vec, limits, 'Color', [0,0,0.545]);
                    if i == 1
                        vcurr.DisplayName= 'v(t, x) = 0';
                    else
                        vcurr.HandleVisibility = 'off';
                    end
                end
            end
            
            
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
            
            xlabel('t')
            ylabel('x_1')
            zlabel('x_2')
            title('Certificate of Disconnectedness', 'FontSize', obj.FS_title)
 

            
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
            
            v_sep = @(t,x,y) obj.out.func.v([t;x;y]);
            v_sep_vec = @(t,x,y)cell2mat(arrayfun(@(i)v_sep(t(i),x(i), y(i)),...
    (1:length(x)),'UniformOutput',false));
            
            
            fimplicit3(@(t,x,y) v_sep_vec(t,x,y), limits, 'MeshDensity', 120)
            
%             
       end
       
    end
end

