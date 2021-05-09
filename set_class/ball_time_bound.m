function [Tmax] = ball_time_bound(n, d)
%BALL_TIME_BOUND returns the maximum length of the geodesic curve present
%in a semialgebraic set {Unit ball} intersect {f(x) >= 0} where x is a
%polynomial in n variables with degree d. This is the supremal length over
%multiple connected components

%bound from Theorem 10.6 of []

%is this an upper bound or a lower bound? important
%the paper makes it seem like an upper bound, which is insufficient for our
%needs
Tmax = (2*d)^(n-1) * n^(-0.5);
end

