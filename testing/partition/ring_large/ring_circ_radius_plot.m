% load('ring_small_2_func.mat')

figure(77)
clf
v_sep = @(t,r) out.func.v([t;r]);
v_sep_vec = @(t,r)cell2mat(arrayfun(@(i)v_sep(t(i), r(i)),...
    (1:length(r)),'UniformOutput',false));

limits =  [0, opt.Tmax, 0, 1];
fsurf(@(t,r) v_sep_vec(t,r), limits,'DisplayName','v(t, x)')
    
ylabel('radius')
xlabel('time')
zlabel('v')

hold on
scatter3(0,X0,out.func.v([0, X0]), 400, 'ko', 'DisplayName', 'X0', 'LineWidth', 3)
scatter3(opt.Tmax,X1,out.func.v([opt.Tmax, X1]), 400, 'k*', 'DisplayName', 'X1', 'LineWidth', 3)

title('Auxiliary Function on Ring-Circle (radius=0.6)', 'FontSize', 16)
legend('location', 'northwest')
% fsurf(@(t,x) func.v(t, x),)