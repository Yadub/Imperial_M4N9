function [ Lambda ] = QR_rayleighshift( A , tol)
%QRSHIFT Computes the eigenvalues and eigenvectors of a triadiagonal matrix
%   Input: A: NxN tridiagonal matrix
%          tol: level of tolerance accepted when implementing shifts 
%   Output: Lambda: Nx1 vector of approximated eigenvalues of A

[N,~] = size(A);    % Set N
Lambda = zeros(N,1);% Initialize vector to store eigenvalues

for k = N:-1:2
    [n,~] = size(A); % Set size of submatrix of A
    while A(k-1,k) > tol
        [Q,R] = mgs_t(A - A(k,k)*eye(n)); % Compute the QR factorization
        A = R * Q + A(k,k) * eye(n); % Iterate A
%         RQ = RQ_givens(A - A(k,k)*eye(n)); % Compute the QR factorization
%         A = RQ + A(k,k) * eye(n); % Iterate A
    end
    Lambda(k) = A(k,k); % Store eigenvalue
    A = A(1:k-1,1:k-1); % Reduce the dimension of A to solve the rest of the problem
end

Lambda(1) = A(1,1); % Store the final eigenvalue

end

