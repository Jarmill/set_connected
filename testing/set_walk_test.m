
s_opt = set_sample_options;
% x0 = @() [0;0];

s_opt.x0 = @()ball_sample(1,2)';

s_opt.Tmax = 4;
s_opt.dt = 0.1;
s_opt.X_func = @blank_event;

rng(33, 'twister')

Np = 10;

out_sim = set_walk(Np, s_opt);

figure(1)
clf
hold on 
for i = 1:Np
%     out_sim = set_walk(x0(), X_func, @() u_func(2), Tmax, dt);
plot(out_sim{i}.x(1, :), out_sim{i}.x(2, :), 'k');
end

function ucurr = u_func(n)
    uraw = 2*rand(n, 1)-1;
    ucurr = uraw/max(abs(uraw));
end
