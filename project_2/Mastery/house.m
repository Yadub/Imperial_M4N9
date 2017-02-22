function [W, A] = house(A)

% This function performs the QR decomposition using the Householder reflections.  
% The input A is overwritten to give the matrix R.
% Input: A: complex mxn matrix with m>=n
% Outputs: W: mxn lower-triangular matrix containing the vectors v defining
%             the Householder reflections.
%          A: nxn upper-triangular matrix containing the coefficients to write
%             columns of A as a linear combination of the columns of Q (to be determined from W).

[m, n] = size(A);
W = zeros(m,n); % preallocate memory for W

for k = 1:n
    V = zeros(m - k + 1,1);
    V = A(k:m, k); % set the vector x (stored as V)
    
    nV = norm(V); %norm of V
    V1 = V(1); %store the first entry of V
    
    if(V1 == 0)
        V(1) = nV + V1; 
    else
        V(1) = sign(V1)*nV + V1; % modify the first entry of V
    end
    
    nV = norm(V); %find the norm of V
    V = V/nV; %normalize V to obtain a unit vector.
    
    A(k:m,k:n) = A(k:m,k:n) - (2.0)*V*(V'*A(k:m,k:n)); %reflect the relevant portions of the columns of A
    W(k:m,k) = V; %store V as the kth column of W
end 