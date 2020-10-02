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
        
        
        %use time independent formulation
        %by default, will find some controller (no objective)
        %set_path_feas_box_signed will find an L1-optimal controller
        time_indep(1,1) logical = false; 
        
        
        %objective to minimize:
        %   time:       time to reach a point on X1 from X0
        %   int:        Integral of ||u||^2_2 along trajectories (on box)
        %               or
        objective = 'time';
        
        %function dynamics
        %f: dynamics
        %X: space on which dynamics are valid (arbitrary switching)
%         dynamics = struct('f', [], 'X', {}, 'Tmin', [], 'Tmax', [], 'discrete', 0)
        
        
        %objective to minimize
        %could be a function (single objective)
        %or a cell of functions (minimum of entries)        
%         obj = [];
        
        %iterative cuts
        %output of prior problem may be infeasible (some matrices not PSD)
        %find a feasible point with a cost of approximately prev_cost
%         prev_cost = [];
        
        %% Variables and descriptors
        %variables defining sets (array of symbolic variables)
        x = [];
        
        
        
        
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
        U = 'circ' 
        
        %Coordinate ranges for variables for scaling
        %state variables lie in a box (utils/box_process)
        
        %box is scalar B:           -B      <= x_i <= B
        %box is [Bmin, Bmax]:       -Bmin   <= x_i <= Bmax
        %box is [B_i]:              -Bi     <= x_i <= B_i
        %box is [Bmin_i, Bmax_i]    -Bmin_i <= x_i <= Bmax_i
        box = 1; 
        
        scale = 0; %should variables be scaled to [-1,1] (state) and [0,1] (time)
                        
%         param = [];  %parameters w  
        
        %% additional options
        solver = 'mosek';
        
        %what is the tolerance for identifying a matrix as rank-1?
        rank_tol(1,1) double{mustBePositive} = 1e-3; 
        
        
    end
end