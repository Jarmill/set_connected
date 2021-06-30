% load('ring_small_2_func.mat')

figure(78)
clf
hold on
% v_sep = @(t,r) out.func.v([t;r]);
% v_sep_vec = @(t,r)cell2mat(arrayfun(@(i)v_sep(t(i), r(i)),...
%     (1:length(r)),'UniformOutput',false));

% limits =  [0, 2, 0, 1];
% fsurf(@(t,r) v_sep_vec(t,r), limits,'DisplayName','v(t, x)')

fcell = out.func.func_cell;
for i = 1:length(fcell)
    curr_v = fcell{i}.v;
    curr_lim = reshape(fcell{i}.box', 1, []);
    fsurf(@(t,r) curr_v([t/opt.Tmax;r]), curr_lim,'DisplayName','v(t, x)')
    fcontour(@(t,r) curr_v([t/opt.Tmax;r]), curr_lim, 'k', 'LineWidth', 4, 'DisplayName','v(t, x)=0', 'LevelList', 0)
end

ylabel('radius')
xlabel('time')
zlabel('v')

hold on
scatter3(0,X0,out.func.v([0, X0]), 400, 'ko', 'DisplayName', 'X0', 'LineWidth', 3)
scatter3(2,X1,out.func.v([2, X1]), 400, 'k*', 'DisplayName', 'X1', 'LineWidth', 3)

title('Auxiliary Function on ring-circle')
% fsurf(@(t,x) func.v(t, x),)