function [b] = algo_p1(A,b)
% Yadu Bhageria
% 00733164

% Solves A*u = b using householder triangularization and enforcing the zero-mean condition
% Input: A: mxn matrix specified in the question
%        b: mx1 column vector containing values of f(x) with a periodic f
% Outputs: b: mx1 complex column vector containing the solution

[m ,n] = size(A); % Compute the dimensions of the input

[W, R] = house(A); % Compute householder matrices W and R

% Calculation of Q* b
for k = 1:n
    b(k:m) = b(k:m) - W(k:n,k) * 2.0 * ( W(k:n,k)' * b(k:m) );
end

% Caclulation of R* b i.e. backwards substitution
for i=n:-1:1
   b(i) = ( b(i) - R(i,i+1:n) * b(i+1:n) ) / R(i,i);
end

end