function [Z] = construct_Z(x,N)
% Yadu Bhageria
% 00733164

% Takes an input vector x and interger N that is the size of x. Teturns a matrix Z with the eigenvectors of D2

sqN = sqrt(N); % Saves computing this for each iteration

Z = zeros(N); % Initalize matrix
for k = 0:N-1
    Z(:,k+1) = exp( 1i * ( k ) * x) / sqN; % Compute columns
end

end