classdef coef_con
    %COEF_CON Store coefficients and constraints for set connectedness
    %   Detailed explanation goes here
    
    properties
        coef;
        con;
    end
    
    methods
        function obj = coef_con(coef,con)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.coef = coef;
            obj.con = con;
        end
        
        function obj = append(obj,cc)
            %Append another cc to the current cc (coef_con)
            %   Detailed explanation goes here
            obj.coef = [obj.coef; cc.coef];
            obj.con =  [obj.con; cc.con];
        end
    end
end

