function y = left(x)
% @SUPCON/LEFT - Extract left hand-side

% D. Henrion, January 21, 2004
% Last modified on 16 March 2006

[nr,nc] = size(x);
y = mpol(zeros(nr,nc));
for r = 1:nr
 for c = 1:nc
  y(r,c) = x(r,c).left;
 end
end

