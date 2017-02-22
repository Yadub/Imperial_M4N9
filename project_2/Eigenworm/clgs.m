function [Q, R] = clgs(A)

% This function performs the reduced QR decomposition using the classical
% Gram-Schmidt algorithm.
% Input: A: complex mxn matrix with m>=n
% Outputs: Q: mxn matrix containing the orthonormal vectors spaning
%             range(A)
%          R: nxn upper-triangular matrix containing the coefficients to write
%             columns of A as a linear combination of the columns of Q.

[m, n] = size(A);

for j = 1:n
    V = A(:,j); % initialize vector V
    AV = V;     % initialize vector AV as the jth column of A
    
    for i = 1:j-1
        q = Q(:,i); % set vector q as the ith column of Q
        rij = q'*AV; % take inner product of q and AV
        V = V - q*rij; % subtract component of AV in direction q from V.
        R(i,j) = rij; % store value of the inner product in R
    end
    
    rjj = norm(V); % compute the length of V
    R(j,j) = rjj; % store length as the diagonal of R
    Q(:,j) = V/rjj; % compute and store the jth column of Q from V.
end