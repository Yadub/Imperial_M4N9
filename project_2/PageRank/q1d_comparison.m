% Yadu Bhageria
% CID: 00733164

format shorte;  % Set format
load('G.mat');  % Load given data
alpha = 0.85;   % Set alpha
tol = 1e-8;     % Set tolerance

[N,~] = size(G);% Determine N

x0 = zeros(N,1);% Initalize x0
x0(1) = 1;      % Set x0 = e1

alphas = 0:0.05:0.95;    % Initialize array of alpha values to iterate over
na = length(alphas);    % # alphas to iterate over
tt = zeros(1,na);       % Initialize vector to store the iteration count
tt_sparse = zeros(1,na);% Initialize vector to store the eigenvalues, lambda

for i = 1:na
    alpha = alphas(i);
    
    tic;
    PageRank(G,alpha,tol,x0); % Time the PageRank algorithm
    tt(i) = toc;
    
    tic;
    sparsePageRank(G,alpha,tol,x0); % Time the sparse PageRank algorithm
    tt_sparse(i) = toc;
end

mainfig = figure;
plot(alphas,tt,'-x',alphas,tt_sparse,'-x');
legend('Standard Method','Sparse Method');
title('Comparing time taken by both methods over varying alpha');
xlabel('alpha');
ylabel('Time Taken');
