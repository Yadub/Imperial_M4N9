% Yadu Bhageria
% CID: 00733164

m = 128; % Number of x points
n = 8; % Number of polynomials to be computed

% Construct equally spaced points x_i
delta_x = 2/m;
x = zeros(m,1);
for i = 1:m
    x(i) = delta_x / 2 + (i - 1) * delta_x - 1;
end

Q = chebyshev_mmgs(x, n);

% Plot results
clf;
hold on;
for i = 1:n
    plot( x, Q(:,i))
end
xlabel('x');
title(['Yadu Bhageria: Project 0 part a. n = ' num2str(n) ', m = ' num2str(m)]);
hold off;