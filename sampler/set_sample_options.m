classdef set_sample_options
    %SET_SAMPLE_OPTIONS options for the set_walk sampler of set
    %disconnectedness.
    %   Detailed explanation goes here
    
    properties
        Tmax = 1;   %terminal time
        dt = 0.05;  %increment time in random walk 
        
        %ode function handle
        odefcn = @ode15s;
        %support handle
        X_func = @blank_event;
        
        %nonnegative function handle
        nonneg_func = [];
        
        %function to choose initial point
        x0;
        
        %parallel not yet implemented
        parallel = 0;
    end
    
end

