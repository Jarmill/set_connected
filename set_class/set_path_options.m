classdef set_path_options < handle
    %attributes of peak set_path_options routine
    %a default set of options
    %
    %Does there exist a path from X0 to X1 entirely within X?
    %
    properties
        
        %% properties of run
        %terminal time   
        Tmax(1,1) double{mustBePositive}  = 5;           
        
        %threshold for initial set, v >= epsilon
        epsilon(1,1) double{mustBePositive}  = 1;           
        
        
        %% Variables and descriptors
        %variables defining sets (array of symbolic variables)
        
        t = []; %time
        x = []; %state
        
        verbose = 0; %solver output: https://yalmip.github.io/faq/runinsilent/
        
        %% support sets
        %type @mpol/supcon
        %(X, T): Xt is the support of trajectories at time Tt
%         state_init = []; %initial state at time 0
%         state_supp = []; %states to consider
%         R = 10; %sum(x.^2) <= R^2
        X = [];     %valid set
        X0 = [];    %initial set
        X1 = [];    %terminal set
        
        %options for controller u: 
        %   'circ':     u'u <= 1
        %   'box':      -1 <= ui <= 1
        %   'box_aff':  'box' with new measures
%         U = 'circ' 
        
        %Coordinate ranges for variables for scaling
        %state variables lie in a box (utils/box_process)
        
        %box is scalar B:           -B      <= x_i <= B
        %box is [Bmin, Bmax]:       -Bmin   <= x_i <= Bmax
        %box is [B_i]:              -Bi     <= x_i <= B_i
        %box is [Bmin_i, Bmax_i]    -Bmin_i <= x_i <= Bmax_i
        box = 1; 
        
        scale = 1; %should time be scaled to [0,1] 
                        
%         param = [];  %parameters w  
        
        %% additional options
        solver = 'mosek';
        
        %what is the tolerance for identifying a matrix as rank-1?
%         rank_tol(1,1) double{mustBePositive} = 1e-3; 
        
        
    end
end