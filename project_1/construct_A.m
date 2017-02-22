function [A] = construct_A(N)
% Yadu Bhageria
% 00733164

% Takes an input integer N and returns a matrix A as specified in project 1

dx = 2*pi/ N; % Compute dx
dx2 = dx*dx; % Saves computing this for each iteration

A = zeros(N); % Initalize matrix
% Compute values as specified
for i = 2:N
   A(i,i) = 2;
   A(i,i-1) = -1;
   A(i-1,i) = -1;
end
A(N,1) = -1;
A(1,:) = 1;
A = A / dx2; 

end