% Yadu Bhageria
% CID: 00733164

format shorte;  % Set format
load('G.mat');  % Load given data
tol = 1e-8;     % Set tolerance

[N,~] = size(G); % Determine N

x0 = zeros(N,1); % Initalize x0
x0(1) = 1; % Set x0 = e1

alphas = 0:0.05:0.95;  % Initialize array of alpha values to iterate over
na = length(alphas);    % # alphas to iterate over
X = zeros(N,na);        % Initialize matrix to store outputted x states
C = zeros(1,na);        % Initialize vector to store the iteration count
L = zeros(1,na);        % Initialize vector to store the eigenvalues, lambda
S = zeros(1,na);

for i = 1:na
    alpha = alphas(i);
    [x, count] = PageRank(G,alpha,tol,x0); % Perform the PageRank algorithm
    X(:,i) = x;         % Store outputted state x
    C(i) = count;       % Store iteration count
    P = construct_P(G,alpha); % Construct P for computing the eigenvalues
    S(i) = (x' * x);
    L(i) = x' * P' * x / S(i); % Compute and store the corresponding eigenvalue
end

mainfig = figure;

subplot(2,2,1);
plot(alphas,C,'-x'); % To visualize number if iterations needed for each value of alpha
xlabel('alpha value');
ylabel('# Iterations');
title('# Iterations needed for convergence to 1e-8 for varying alpha');

subplot(2,2,2);
maxvals = max(X,[],1);
plot(alphas,maxvals,'-x');
xlabel('alpha value');
ylabel('Maxval(x)');
title('Maxval(x) for varying alpha');

subplot(2,2,3);
plot(alphas,L,'-x');
xlabel('alpha value');
ylabel('Eigenvalue');
title('Eigenvalue for the outputted x state for varying alpha');
axis([0 1 1-1e-8 1+1e-8]);

subplot(2,2,4);
plot(alphas,S,'-x');
xlabel('alpha value');
ylabel('Size of x^T x');
title('Size of x^T x for varying alpha');