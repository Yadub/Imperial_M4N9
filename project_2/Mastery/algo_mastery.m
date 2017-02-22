function [ x ] = algo_mastery( A, b )
%ALGO_MASTERY Summary of this function goes here
%   Detailed explanation goes here

% Solves (P^T D2 P) (P^T x) = P^T b using Band Gaussian elimination
% Input: D2: NxN matrix specified in the question
%        P: NxN permutation matrix specified in the question
%        b: Nx1 column vector containing values of f(x) with a periodic f
% Outputs: x: Nx1 column vector containing the solution

[N,~] = size(A); % Set N

p = 2; % Lower band size
q = 2; % Upper band size

% Use LU decomposition Algo
for k = 1:N-1
    for i = k+1:min(k+p,N)
       A(i,k) = A(i,k)/A(k,k);
    end
    for j = k+1:min(k+q,N)
        for i = k+1:min(k+p,N)
            A(i,j) = A(i,j) - A(i,k)*A(k,j);
        end
    end
end

% Forward substitution i.e. L' b0. Note L' has 1's on the middle diagonal
b(2) = b(2) - A(2,1) * b(1);
for k = 3:N
    b(k) = b(k) - A(k,k-2:k-1) * b(k-2:k-1);
end

% Caclulation of U' b i.e. backwards substitution. Note U(N,N) is assumed
% to be 1
b(N-1) = ( b(N-1) - A(N-1,N) * b(N) ) / A(N-1,N-1);
for i=N-2:-1:1
   b(i) = ( b(i) - A(i,i+1:i+2) * b(i+1:i+2) ) / A(i,i);
end

x = b;

end

