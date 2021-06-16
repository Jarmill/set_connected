load('ring_circle_large_2_func.mat')
Tmax = 2;
X0 = 0.8;
X1 = 0.3;

figure(77)
clf
v_sep = @(t,r) func.v([t;r]);
v_sep_vec = @(t,r)cell2mat(arrayfun(@(i)v_sep(t(i), r(i)),...
    (1:length(r)),'UniformOutput',false));

limits =  [0, Tmax, 0, 1];
fsurf(@(t,r) v_sep_vec(t,r), limits,'DisplayName','v(t, x)')
    
ylabel('radius')
xlabel('time')
zlabel('v')

hold on
scatter3(0,X0,func.v([0, X0]), 400, 'ko', 'DisplayName', 'X0', 'LineWidth', 3)
scatter3(Tmax,X1,func.v([Tmax, X1]), 400, 'k*', 'DisplayName', 'X1', 'LineWidth', 3)

title('Auxiliary Function on Ring-Circle (radius=0.6)', 'FontSize', 16)
legend('location', 'northwest')
% fsurf(@(t,x) func.v(t, x),)


%% level set
figure(78)
clf
hold on
scatter(0,X0, 400, 'ko', 'DisplayName', 'X0', 'LineWidth', 3)
scatter(Tmax,X1, 400, 'k*', 'DisplayName', 'X1', 'LineWidth', 3)

for i = 1:length(out_sim)
    osc = out_sim{i};
    r_curr = sqrt(sum(osc.x_traj.^2, 2));
    if i==1
        plot(osc.t_traj, r_curr, 'c', 'DisplayName', 'Trajectory')
    else
        plot(osc.t_traj, r_curr, 'c', 'HandleVisibility', 'off')
    end
end
fimplicit(@(t,r) v_sep_vec(t,r), limits,'DisplayName','v(t, x)=0', 'Color', 'k', 'LineWidth', 3)

plot([0,2], [0.6, 0.6], '--', 'color', 0.6*[1,1,1], 'LineWidth', 3, 'DisplayName', 'Boundary')
plot([0,2], [0.7, 0.7], '--', 'color', 0.6*[1,1,1], 'LineWidth', 3, 'HandleVisibility', 'off')
ylim([0,1])
xlim([0,Tmax])
legend('location', 'southwest')
title('Certificate Contour', 'fontsize', 14)
% end

