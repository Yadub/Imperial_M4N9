function [V] = chebyshev_mmgs(x, n)
% Computes N chebyshev polynomials over m equally spaced x points

[m, ~] = size(x);
% Construct the resulting Vandermonde matrix
V = zeros(m,n); 
for i = 1:n
    V(:,i) = x.^(i-1);
end
% Compute the first n Chebyshev polynomials and store them in V
R = zeros(m,n); % Initialize R
w = (1 - x.^2).^(-1/2); % Compute the values of the weight function
for i = 1:n
    R(i,i) = ( V(:,i)' * ( V(:,i) .* w)) ^ (1/2); % modified norm 
    % It is worth noting that I do not multiply by delta_x in the norm as we later normalize the polynomials.
    q = V(:,i)/R(i,i);
    for j = i+1:n
        R(i,j) = q' * ( V(:,j) .* w); % modified norm
        V(:,j) = V(:,j) - R(i,j) * q;
    end
    V(:,i) = q;
end
%Find the values of the polynomials at the last x value
T = zeros(n,1); 
T(1) = 1;
T(2) = x(m);
for i = 3:n
    T(i) = 2 * x(m) * T(i-1) - T(i-2);
end
% Normalize the values of Q
for i = 1:n
    scaling_factor = T(i)/V(m,i);
    V(:,i) = V(:,i) * scaling_factor;
end
