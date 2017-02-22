function [ lambda ] = construct_lambda( N )
% Yadu Bhageria
% 00733164

% Returns an array of size N with the eigenvalues of D_2
%   Input: Integer N with the number of eigenvalues wanted
%   Output: Array containing the N eigenvalues of D_2

dx = 2*pi/N; % Compute dx
dx2 = dx*dx; % Saves computing this for each iteration

lambda = zeros(N,1); % Initalize array
for k = 2:N
        lambda(k) = 2 * ( cos( (k-1) * dx ) - 1) / (dx2); % Compute values
end

end

