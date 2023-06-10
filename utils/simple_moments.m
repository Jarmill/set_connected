function [mom_out] = simple_moments(var_count, order)
%SIMPLE_MOMENTS monomial powers and moments without symmetry
%Inputs:
%  var_count: an array of length 2
%               [symmetric;
%                standard]
%
%           the symmetry is f(xs, x) = f(-xs, x) 
%           (symmetric, standard)
%
% order:    degree of relaxation, monomial powers
%           moment matrix will have moments of up to degree d=2*order
%
%Outputs:   struct mom_out with fields:
%   monom_int:      monom_int: monomials that are on the interior of the
%                   moment matrices
%   index:          Don't know what this does
%   M:              indices for moments (vectorized, not a matrix)
%   N:              Number of moments on the interior of the moment matrix
%   order:          order of relaxation
%% preliminary counting and indexing

d = 2*order;

mom_out = struct;

%% Form monomial basis 

%find the monomial basis
[basis, ~] = momentPowers(0, sum(var_count), order);

%split monomial basis into even and odd by symmetric variables
% basis_sym = sum(basis(:, ind_sym), 2);
% basis_mod = mod(basis_sym, 2);
% 
% basis_even = basis(basis_mod==0, :);
% basis_odd  = basis(basis_mod==1, :);


%% Find monomial powers inside the matrix
% m_even = size(basis_even, 1);
% temp_even = bsxfun(@plus,kron(ones(m_even,1),basis_even),kron(basis_even,ones(m_even,1)));
% 
% 
% m_odd= size(basis_odd, 1);
% temp_odd = bsxfun(@plus,kron(ones(m_odd,1),basis_odd),kron(basis_odd,ones(m_odd,1)));
m = size(basis, 1);
temp = bsxfun(@plus,kron(ones(m, 1), basis), kron(basis, ones(m, 1)));

%monomials indexed by moment matrix
monom_int = unique(temp, 'rows');
index = 1:size(monom_int, 1);

% index_order2 = any(monom_int(:, 1:sum(var_count(1:2))) == 2, 2);

%vectorized moment matrix
M = get_vec_ind_func(monom_int,temp,index,0)';



% ZZ_even = unique(temp_even, 'rows');
mom_out.monom_int = monom_int;
mom_out.index = index;
% mom_out.index_order2 = index_order2; 
% mom_out.M_even = M_even;
mom_out.M = M;

%number of monomials
mom_out.N_int = index(end);
mom_out.N= m;
% mom_out.N_odd = m_odd;
mom_out.order = order;

end

