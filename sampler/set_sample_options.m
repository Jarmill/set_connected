classdef set_sample_options
    %SET_SAMPLE_OPTIONS options for the set_walk sampler of set
    %disconnectedness.
    %   Detailed explanation goes here
    
    properties
        Tmax = 1;
        dt = 0.05;
        
        
        odefcn = @ode15s;
        
        X_func = @blank_event;
        
        x0;
        
        parallel = 0;
    end
    
end

