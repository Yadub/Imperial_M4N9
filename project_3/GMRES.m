function [ b, iters] = GMRES( A, b, tol, preconditioning )
%GMRES Summary of this function goes here
%   Detailed explanation goes here

[N,~] = size(A); % Set N

% PRECONDITIONING
if preconditioning
    [L,U] = ilu(A); % Compute incomplete LU factorization
    b = U \ ( L \ b);    % Use it to precondition starting b
end

H = zeros( N+1, N);    % Initialize H
Q = zeros( N, N+1);    % Initialize Q
s = zeros(N,1);        % Initialize s vector for givens rotations
c = zeros(N,1);        % Initialize c vector for givens rotations

nb = norm(b);   % Compute the norm of b
Q(:,1) = b / nb;% Set b as the first normalized basis vector
b(1) = nb;      % Vector b is used to store y from now on

for i = 1:N
    v = A * Q(:,i); % Sparse matrix vector product
    
    if preconditioning
        v = U \ ( L \ v);
    end
    
    % Arnoldi Iteration
    for j = 1:i
        H(j,i) = Q(:,j)' * v;   
        v = v - H(j,i) * Q(:,j);% MGS
    end
    H( i+1, i) = norm(v);       
    Q( :, i+1) = v / H( i+1, i);
    
    % Perform Given's Rotations on new values
    for j = 1:i-1
        Hj = H(j,i);        % Store temp var
        Hjp1 = H( j+1, i);  % Store temp var
        H(j,i) = c(j) * Hj - s(j) * Hjp1;
        H( j+1, i) = s(j) * Hj + c(j) * Hjp1;
    end
    % Compute coefficients for new Givens Rotation that sets H(i+1,i) to 0
    r = sqrt( H(i,i)^2 + H( i+1, i)^2);
    c(i) = H(i,i) / r;
    s(i) = - H( i+1, i) / r;
    % Compute values using the rotation
    H(i,i) = r;             % r = c(i)*H(i,i) - s(i)*H(i+1,i)
    H( i+1, i) = 0;         % By definition of the rotation being performed
    b(i+1) = s(i) * b(i);   % Using the fact that b(i+1) = 0
    b(i) = c(i) * b(i);     % Using the fact that b(i+1) = 0
    % Break loop if relative tolerance is achieved
    if abs( b(i+1) / nb) < tol
        break;
    end    
end

% Backwards substitution
for j = i:-1:1
    b(j) = b(j) - H(j,j+1:i) * b(j+1:i);
    b(j) = b(j) / H(j,j);
end
  
b = Q(:,1:i) * b(1:i);  % Compute solution
iters = i;              % Set number of iterations to be outputted
end

