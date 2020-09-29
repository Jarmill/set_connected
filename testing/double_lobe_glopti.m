%try to do proofs of disconnectedness of double-lobe region
%with gloptipoly

%Points in left lobe
%[-1; 0]
%[-0.75; -0.25]
%[-1.25; 0.25]

% x0 = [-0.75; -0.25];
% x1 = [-1.25; 0.25];
% x1 = [-1; 0];

%Points on right lobe
% x0 = [1.25; -1];
%[1; 1]
% x1 = [1.5; 0.5];


FEAS = 1;

if FEAS
    x0 = [1.25; -1];
    x1 = [1.5; 0.5];
else
    x0 = [-0.75; -0.25];
    x1 = [1.5; 0.5];
end



%constraint set
f = @(x) -(x(1)^4 + x(2)^4 - 3*x(1)^2 - x(1)*x(2)^2 - x(2) + 1);


n = 2;

order = 1;
d =2*order;
T = 2; %maximum time


%% Start the variables and measures
% mpol('t00', 1, 1);
% mpol('x00', n, 1);
% mu0 = meas([t00; x00]);

mpol('tT', 1, 1);
mpol('xT', n, 1);
muT = meas([tT; xT]);

mpol('t', 1, 1);
mpol('x', n, 1);
mpol('u', n, 1);
mu = meas([t; x; u]);

%support constraints
% supp_mu0 = [t00 == 0; x00 == x0];
supp_muT = [tT*(T - tT) >= 0; xT == x1];
supp_mu = [t*(T - t) >= 0; f(x) >= 0; 1 - u'*u >= 0];

%moment constraints
% mon0 = mmon([t00; x00], d);
monT = mmon([tT; xT], d);
mon  = mmon([t; x], d);

yT = mom(monT);
Ay = mom(diff(mon, x)*u) + mom(diff(mon, t));
       

powers = genPowGlopti(n+1, d);
y0 = prod(([0; x0]').^powers, 2);



Liou = yT - Ay - y0;

%altogether now
%supp_con = [supp_mu0; supp_muT; supp_mu];
supp_con = [supp_muT; supp_mu];
mom_con  = (Liou == 0);

objective = min(mom(tT));

mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));

P = msdp(objective, ...
    mom_con, supp_con);

[status,obj_rec, m,dual_rec]= msol(P);

%Read output
% M0 = double(mmat(mu0));
dual_rec_v = dual_rec{1};
v = dual_rec_v'*mon;
v1 = subs(v, x, x1);

vval = @(ti, xi) eval(v, [t; x], [ti; xi]);    %dual v(t,x,w)
Lvval = @(ti, xi) eval(Lv, [t; x], [ti; xi]);   %Lie derivative Lv(t,x,w)
v1val = @(ti) eval(v1, [t], [ti]);    %dual v(t,x,w)
       

% v0 = eval(v, [t; x], [0; x0]);
% v1t = subs(v, x, x1);

if status == 0
    %feasible problem
    powers_half = genPowGlopti(n+1, order);
    y0_half = prod(([0; x0]').^powers_half, 2);
    M0 = y0_half*y0_half';


    MT = double(mmat(muT));
    Mocc = double(mmat(mu));
elseif status == -1
    %infeasible problem
end
tlist = linspace(0, T, 200);
tvt = v1val(tlist);
plot(tlist, tvt+tlist)
v0 = vval(0, x0);

