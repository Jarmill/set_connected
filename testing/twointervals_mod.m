% Liouville certificate of disconnectedness of the union of two intervals
% using YALMIP SOS modeling

% D. Henrion, 22 Apr 22
d = 2; % SOS degree

sdpvar t x u

e = 1/4; 
e=0.01;
g = -(1+x)*(-e-x)*(-e+x)*(1-x); % X=[-1,-e]U[e,1]
g0 = -(x+(1+e)/2)^2; % X0={-(1+e)/2}
g1 = -(x-(1+e)/2)^2; % X1={(1+e)/2}

[v,vc,vb] = polynomial([t,x],d);
[st,stc] = polynomial([t,x,u],d-2);
[sx,sxc] = polynomial([t,x,u],d-2);
[su,suc] = polynomial([t,x,u],d-2);
[s0,s0c] = polynomial(x,d-2);
[s1,s1c] = polynomial(x,d-2);

dv = -jacobian(v,t)-jacobian(v,x)*u; % derivative
[dvc,dvb] = coefficients(dv,[t,x,u]);

solvesos([sos(dv-st*t*(1-t)-sx*g-su*(1-u^2)); ... % -dv > 0 everywhere
    sos(-replace(v,t,0)-1-s0*g0); ... % -v-1 > 0 at t=0 and X0
    sos(replace(v,t,1)-s1*g1); ... % v > 0 at t=1 and X1
    sos(st); sos(sx); sos(su); sos(s0); sos(s1)],...
    [norm(vc)],[],[vc;stc;sxc;suc;s0c;s1c]);

vs = vectorize(sdisplay(double(vc)'*vb));
dvs = vectorize(sdisplay(double(dvc)'*dvb));
[t,x] = meshgrid([0,1],[-5/4,5/4],1000);

% plot (t,x) -> v(t,x)
close all
surf(t,x,eval(vs));
xlabel t
ylabel x
zlabel v

figure
% plot (t,x) -> dv(t,x,u) for different values of u
h = []; hl = {};
for u = linspace(-1,1,10)
 h = [h; surf(t,x,eval(dvs))];
 hl = {hl{:},['u=' num2str(u)]};
 hold on;
end
xlabel t
ylabel x
zlabel dv
legend(h,hl)
