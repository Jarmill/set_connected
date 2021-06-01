classdef set_manager_partition < handle
    %SET_MANAGER_PARTITION 
    %Given sets X0 and X1 within a base space X,
    %determine if there exists a path between points in X0 and X1 (feas)
    %or there does not exist such a path (infeas)
    properties
        options;
        spacing; %number of locations at each dimension
        loc;     %the locations on the dimension
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
                obj.spacing = spacing;
            end
            
            obj.loc = cell(obj.spacing);
            
            if isempty(obj.options.t)
                obj.options.t = sdpvar(1, 1);
            end
        end
        s
    end
end

