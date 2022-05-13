% Liouville certificate of disconnectedness of the union of two intervals
% using YALMIP SOS modeling

% D. Henrion, 22 Apr 22

d = 2; % SOS degree
% d=6;

sdpvar t x u


% e = 1/4; 
% e = 0.1;
e = 0.05;
% e = 0.025;
% e=0.01;
epsilon = 1;

% xbox = [-1;1];      %works, feasible barrier
xbox = [-0.1;1];     %doesn't work, infeasible SDP


T = 1;
offset = 0;
% low =  offset-e;
% high = offset+e;

low_alpha = 0.7;
high_alpha = 0.75;
low = xbox(2)*low_alpha+ xbox(1)*(1-low_alpha);
high = xbox(2)*high_alpha+ xbox(1)*(1-high_alpha);


% X0_alpha = 0.2;
X0_alpha = 0.05;
% X1_alpha = 0.9;
X1_alpha = 0.99;

X0 = xbox(2)*X0_alpha + xbox(1)*(1-X0_alpha);
X1 = xbox(2)*X1_alpha + xbox(1)*(1-X1_alpha);

% X0 = (xbox(1)+low)/2;
% X1 = (xbox(2)+high)/2;


% g = -(1+x)*(-e-x)*(-e+x)*(1-x); % X=[-1,-e]U[e,1]
g = -(x-xbox(1)).*(xbox(2)-x).*(x-low).*(high-x);
gf = polyval_func(g, x);
g0 = -(x-X0)^2; % X0={-(1+e)/2}
g1 = -(x-X1)^2; % X1={(1+e)/2}

[v,vc,vb] = polynomial([t,x],d);
[st,stc] = polynomial([t,x,u],d-2);
[sx,sxc] = polynomial([t,x,u],d-2);
[su,suc] = polynomial([t,x,u],d-2);
[s0,s0c] = polynomial(x,d-2);
[s1,s1c] = polynomial(x,d-2);

dv = -jacobian(v,t)-jacobian(v,x)*u; % derivative
[dvc,dvb] = coefficients(dv,[t,x,u]);

cons = [sos(dv-st*t*(T-t)-sx*g-su*(1-u^2)); ... % -dv > 0 everywhere
    sos(-replace(v,t,0)-epsilon-s0*g0); ... % -v-1 > 0 at t=0 and X0
    sos(replace(v,t,1)-s1*g1); ... % v > 0 at t=1 and X1
    sos(st); sos(sx); sos(su); sos(s0); sos(s1)];

objective = norm(vc);

solvesos(cons,...
    objective,[],[vc;stc;sxc;suc;s0c;s1c]);

cv_rec = double(vc);
Lcv_rec = double(dvc);

v_rec = cv_rec'* monolist([t; x], d);
    nv_rec = norm(cv_rec);
fv = polyval_func(v_rec, [t; x]);
vs = vectorize(sdisplay(double(vc)'*vb));
dvs = vectorize(sdisplay(double(dvc)'*dvb));
[t,x] = meshgrid([0,T],xbox',1000);

%% plotting of barrier function
% plot (t,x) -> v(t,x)
% close all
figure(1)
clf
% surf(t,x,fv([t,x]));
fsurf(@(t,x) fv([t;x]), [0,T, xbox'])
xlabel t
ylabel x
zlabel v
hold on
zl = zlim;
xl_pattern = [0, 1, 1, 0, 0];
zl_pattern = zl([1, 1, 2, 2, 1]);
patch(xl_pattern, low*ones(5,1), zl_pattern, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'None', 'HandleVisibility', 'Off')
patch(xl_pattern, high*ones(5,1), zl_pattern, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'None', 'DisplayName', 'Region')
    fcontour(@(t,x) fv([t;x]), [0,T,xbox'], 'k', 'LevelList', 0, 'LineWidth', 4, 'DisplayName', 'v(t,x)=0');

% if sol.problem == 0
    
scatter3(0,X0, fv([0; X0]), 400, 'ko', 'DisplayName', 'X0', 'LineWidth', 3)
scatter3(T,X1,fv([T; X1]), 400, 'k*', 'DisplayName', 'X1', 'LineWidth', 3)


%% plot of support function
figure(3)
clf
Nxsample = 300;
xrange = linspace(xbox(1), xbox(2), Nxsample);

gx = gf(xrange);
g0 = gf(X0);
g1 = gf(X1);

hold on
plot(xrange, gx)
scatter(X0, g0, 200, 'ko')
scatter(X1, g1, 200, 'k*')
scatter(low, 0, 200, 'ks')
scatter(high, 0, 200, 'ks')
plot(xbox, [0,0], '--k')
xlim(xbox)
% end
% figure(2)
% clf
% % plot (t,x) -> dv(t,x,u) for different values of u
% h = []; hl = {};
% for u = linspace(-1,1,10)
% %  h = [h; surf(t,x,eval(dvs))];
%  hl = {hl{:},['u=' num2str(u)]};
%  hold on;
% end

% xlabel t
% ylabel x
% zlabel dv
% legend(h,hl)
