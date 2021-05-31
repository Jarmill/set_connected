function out_sim = set_walk_single(x0, X_func, u_func, Tmax, dt)
%set_walk: sample a single random walk starting at the point X0


%work on von mises later, do standard random walk for now
%angles are chosen according to a von mises distribution with parameter
%kappa. the inputs are the directional angles normalized to the boundary of
%the Linf ball (at least one input coordinate will have maximum absolute
%value of 1).
%
%this isn't going to be a completely uniform random walk, but hopefully it
%will be good enough.

n = length(x0);



Nstep = ceil(Tmax/dt);


% dtheta= vmrnd(0, kappa, 1, Nstep);


% theta_curr = vmrnd(0,0);  
xcurr = x0;

% theta = wrap2pi(cumsum(dtheta));
x = zeros(length(x0), Nstep);

odeopts = odeset('Events', X_func);

out_sim=struct;
out_sim.u = zeros(n, Nstep);
out_sim.x = zeros(n, Nstep);
out_sim.x(:, 1) = xcurr;

out_sim.t_traj = [];
out_sim.x_traj = [];


rng('shuffle')

for i = 1:Nstep-1
    
    %applied input inside box [-1,1]^n
 
    ucurr = u_func();

   
    fcurr = @(t,x) ucurr;
    
    [tode, xode] = ode45(@(t,x) ucurr, [0,dt],xcurr, odeopts);
    
    xnext = xode(end, :)';
    
    out_sim.u(:, i) = ucurr;
    out_sim.x(:, i+1) = xnext;
    
    out_sim.t_traj = [out_sim.t_traj; tode + dt*(i-1)];
    out_sim.x_traj = [out_sim.x_traj; xode];
        
    
    
    xcurr = xnext;
    
   if tode(end) ~= dt
       break
   end
        
end


out_sim.u = out_sim.u(:, 1:(i+1));
out_sim.x = out_sim.x(:, 1:(i+1));
% theta = linspace(0, 
% if ~isempty(options.nonneg_func)
%     %evaluate nonnegative functions along controlled trajectories
%     
%     Ntraj = length(out_sim.t_traj);
%     out_sim.nonneg_traj = zeros(2*n+1, N_traj);
%     for i = 1:Ntraj
%         out_sim.nonneg_traj(:, i) = options.nonneg_func(out_sim.x_traj(:, i));
%     end
%         
% 
% 
% end
end

