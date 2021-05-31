function [bound] = T_bound(n, d, r)
%T_BOUND maximum length of geodesic path connecting two points in a 
%semialgebraic set. This set is defined by an n-variate polynomial f(x)>=0
% where f has degree d, contained within the ball with radius r

%based on results from 
% Didier D’Acunto and Krzysztof Kurdyka Bounds for gradient trajectories of definable 
%functions with applications to robotics and semialgebraic geometry. 
% 2003 Sep 23.
if nargin < 3
    r = 1;
end

 
%Remark 4.3 of gradient paper
nu = @(n) 2*gamma(0.5) * gamma((n+1)/2)/gamma(n/2);

%Theorem 7.10 of gradient paper
A = @(n, d) nu(n)*( (3*d-4)^(n-1) + 2*(3*d-3)^(n-2));

%proposition 10.4 of gradient paper
bound = r*A(n, d+2);


end

