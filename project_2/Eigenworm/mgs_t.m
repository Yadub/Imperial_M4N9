function [A, R] = mgs_t(A)

% This function performs the reduced QR decomposition using the modified
% Gram-Schmidt algorithm.  The input A is overwritten to give the matrix Q.
% Input: A: complex NxN symmetric tridiagonal matrix
% Outputs: A: NxN matrix containing the orthonormal vectors spaning range(A)
%          R: NxN upper-triangular matrix containing the coefficients to write
%             columns of A as a linear combination of the columns of Q.

[N, ~] = size(A);
R = zeros(N);

for i = 1:N
    
    rii = norm(A(:,i)); %compute diagonal entry of R from ith column of A
    q = A(:,i)/rii; %compute the ith column of Q
    R(i,i) = rii; %store the entry ii in R
    A(:,i) = q; %store the ith column of Q as the ith column of A
    
    for j = i+1:N
        
        V = A(:,j); %set V as the jth column of A
        rij = q'*V; %take the inner product of q and V
        R(i,j) = rij; %store inner product in R
        V = V - rij*q; %subtract off component V in direction q
        A(:,j) = V; %store V as the jth column of A
        
    end
end