function [ evec ] = evec( A, eval )
%EVEC Computes the eigenvector associated with the given eigenvalue for A
%   Input: A: NxN symmertric matrix
%          eval: an eigenvalue of the matrix A
%   Output: evec: Nx1 eigenvector associated with the given eigenvalue

[N,~] = size(A);% Set N

A = A - eval * eye(N);  % Subtract the eigenvalue from the main diagonal of A

% General LU Decomposition Algorithm from Lecture Notes
% The Pivot Matrix P is ommited as LUx = Pb is being solved for b = 0
U = A;      % Initialize U
L = eye(N); % Initailize L
for k = 1:N-1
    % Select pivot to maximize |U(i,k)|
    i = k;
    for j = k+1:N
        if abs(U(i,k)) < abs(U(j,k))
            i = j;
        end
    end
    % Exchange rows
    U([k, i],k:N) = U([i, k],k:N);
    L([k, i],1:k-1) = L([i, k],1:k-1);
    % Perform the rest of the algo
    for j = k+1:N
        L(j,k) = U(j,k) / U(k,k);
        U(j,k:N) = U(j,k:N) - L(j,k) * U(k,k:N);
    end
end

% Backward substitution
x = zeros(N,1);
x(N) = 1;
for k = N-1:-1:1
    for j = k + 1:N
        x(k) = x(k) - U(k, j) * x(j);
    end
    x(k) = x(k) / U(k,k);
end
evec = x / norm(x); % Normalize output eigenvector 

end
