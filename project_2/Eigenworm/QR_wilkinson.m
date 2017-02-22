function [ Lambda ] = QR_wilkinson( A , tol)
%QRSHIFT Computes the eigenvalues and eigenvectors of a triadiagonal matrix
%   Input: A: NxN tridiagonal matrix
%          tol: level of tolerance accepted when implementing shifts 
%   Output: Lambda: Nx1 vector of approximated eigenvalues of A

[N,~] = size(A);    % Set N
Lambda = zeros(N,1);% Initialize vector to store eigenvalues

for k = N:-1:2
    [n,~] = size(A); % Set size of submatrix of A
    while A(k-1,k) > tol
        B = A(k-1:k,k-1:k); % Set B as the lower 2x2 matrix
        d = ( B(1,1) - B(2,2) ) / 2;
        mu = B(2,2) - ( sign(d) * B(2,1)^2 ) / ( abs(d) + sqrt(d^2 + B(2,1)^2) ); % Set mu
        RQ = RQ_givens( A - mu * eye(n) ); % Compute the QR factorization using modified gram schmidt
        A = RQ + mu * eye(n); % Iterate A
%         [Q,R] = mgs_t( A - mu * eye(n) ); % Compute the QR factorization using modified gram schmidt
%         A = R*Q + mu * eye(n); % Iterate A
    end
    Lambda(k) = A(k,k); % Store eigenvalue
    A = A(1:k-1,1:k-1); % Reduce the dimension of A to solve the rest of the problem
end

Lambda(1) = A(1,1); % Store the final eigenvalue

end


