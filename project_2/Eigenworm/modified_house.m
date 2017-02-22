function [W, A] = modified_house(A)

% This function performs the QR decomposition using the Householder reflections.  
% The input A is overwritten to give the matrix R.
% Input: A: complex mxn matrix with m>=n
% Outputs:  A: mxm upper-triangular matrix containing the coefficients to write
%             columns of A as a linear combination of the columns of Q (to be determined from W).
%           W: mxm lower triangular matrix containing the vectors v
%           defining the householder rreflections

[m, ~] = size(A);

W = zeros(m,m); % preallocate memory for W

for k = 1:m-2

    V = A(k+1:m, k); % set the vector x (stored as V)
    
    nV = norm(V); %norm of V
    V1 = V(1); %store the first entry of V
    
    if(V1 == 0)
        V(1) = nV + V1; 
    else
        V(1) = sign(V1)*nV + V1; % modify the first entry of V
    end
    
    nV = norm(V);   % Find the norm of V
    V = V/nV;       % Normalize V to obtain a unit vector.
    W(k+1:m,k) = V; % Store V in W
    
    % METHOD USING SPARSITY AND SYMMETRY
    A(k+1:m,k) = A(k+1:m,k) - (2.0)*V*(V'*A(k+1:m,k)); %reflect first of the columns of A
    A(k,k+1:m) = A(k+1:m,k); % Use symmtery to copy kth row to the kth column
    
    % Rows 1:k have zeros and thus we can ignore them
    T1 = A(k+1:m,k+1:m) * V; % Precompute T1 to save FLOPs
    T2 = 2*T1 - 2*V*(V'*T1); % Precompute T2 to save FLOPs
    
    % Use symmetry to copy rows k+1:m to the columns k+1:m
    for j = k+1:m
        A(j,j:m) = A(j,j:m) - T2(j-k)*V(j-k:m-k)' - V(j-k)*T2(j-k:m-k)';
        A(j:m, j) = A(j, j:m)';
    end
        
%   METHOD FROM LECTURES
%     A(k+1:m,k:m) = A(k+1:m,k:m) - (2.0)*V*(V'*A(k+1:m,k:m)); %reflect the relevant portions of the columns of A
%     A(1:m,k+1:m) = A(1:m,k+1:m) - (2.0)*(A(1:m,k+1:m)*V)*V';
    
%   MUCH SLOWER METHOD BUT DEMONSTRATES THE MAIN IDEA WELL
%     Vp = zeros(m,1);
%     Vp(k+1:m) = V;
%     P = eye(m) - (2.0)*Vp*Vp';
%     P = sparse(P);
%     A = P*A*P;
end 