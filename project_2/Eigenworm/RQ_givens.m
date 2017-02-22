function [ O ] = RQ_givens( A )
%RQ_GIVENS Computes RQ for a tridiagonal matrix A
%   Input: A: NxN tridiagonal symmetric matrix
%   Output: O: NxN matrix with form RQ from the QR factorization of A 

[N,~] = size(A); % Set N

% Initalize vectors
a = diag(A);
b = diag(A,1);
c = zeros(N,1);
s = zeros(N,1);

t = b(1); % For first iteration
for j = 1:N-1
    r = sqrt( a(j)^2 + t^2 ); 
    c(j) = a(j)/r;                      
    s(j) = t/r;
    a(j) = r;
    t = b(j);
    b(j) = t*c(j) + a(j+1)*s(j);
    a(j+1) = -t*s(j) + a(j+1)*c(j);
    if j ~= N-1
        t = b(j+1);
        b(j+1) = t*c(j);
    end
end

% Compute leading and first off diagonal of the matrix RQ from a, b, c, and s 
for i = 1:N-1
    a(i) = a(i) * c(i) + b(i) * s(i);
    b(i) = a(i+1) * s(i);
    a(i+1) = a(i+1) * c(i);
end

% Create matrix to be outputted
O = diag(a) + diag(b,1) + diag(b,-1);

end
