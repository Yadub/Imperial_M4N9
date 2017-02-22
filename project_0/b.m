% Yadu Bhageria
% CID: 00733164

% Consider the substitution x = cos(theta)

m = 128; % Number of x points
n = 8; % Number of polynomials to be computed

% Construct equally spaced points over theta for x = cos(theta)
x = zeros(m,1);
for i = 1:m
    x(i) = cos( pi * ( 2 * i - 1) / ( 2 * m));
end

Q = chebyshev_smgs(x, n);

% Plot results
clf;
hold on;
for i = 1:n
    plot( x, Q(:,i))
end
xlabel('x');
title(['Yadu Bhageria: Project 0 part b. n = ' num2str(n) ', m = ' num2str(m)]);
hold off;