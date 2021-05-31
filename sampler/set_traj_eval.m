function [out_sim] = set_traj_eval(out, out_sim)
%SET_TRAJ_EVAL evaluate nonnegative trajectories along random-walk
%trajectories
sw = out.poly.scale_weight;

%   Detailed explanation goes here

nn = out.func.nonneg;
v = out.func.v;
n = size(out_sim{1}.x, 1);
for i = 1:length(out_sim)
    tcurr = out_sim{i}.t_traj;
    xcurr = out_sim{i}.x_traj';
    Nt = length(tcurr);
    nonneg_curr = zeros(2*n+1, Nt);
    v_curr = zeros(1, Nt);
    for j = 1:Nt
        nonneg_curr(:, j) = nn([tcurr(j)/sw; xcurr(:,j)]);
        v_curr(:, j) = v([tcurr(j)/sw; xcurr(:,j)]);
    end
    out_sim{i}.v = v_curr;
    out_sim{i}.nonneg = nonneg_curr;
end

end

