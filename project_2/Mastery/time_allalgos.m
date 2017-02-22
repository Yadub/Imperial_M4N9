function [tp1 , tp4, tfft, tm] = time_allalgos(N)
% Yadu Bhageria
% 00733164

% Takes as an input number of points N and returns the time taken by each
% of the three algos implemented

% Construct x points
dx = 2*pi/N;
x = zeros(N,1);
for j = 0:N-1
    x(j+1) = (j) * dx;
end
% Compute f(x)
fx = sin(x);
% Construct the array storing the eigenvalues of D_2
lambda = construct_lambda(N);

if N < 2^11
    
    % Construct the A matrix for solving Au = f(x)
    A = construct_A(N);

    tic;
    algo_p1(A,fx); % Time algo p1
    tp1 = toc;

    % Construct the Z matrix for solving - (Z Lambda Z') u = f(x)
    Z = construct_Z(x,N);

    tic;
    algo_p4(Z,lambda,fx,N); % Time algo p4
    tp4 = toc;
else
    tp1 = 0;
    tp4 = 0;
end
    
tic;
algo_fft(N,lambda,fx); % Time algo fft
tfft = toc;

% Construct D2
D2 = diag(2*ones(N,1)) + diag(-1*ones(N-1,1),1) + diag(-1*ones(N-1,1),-1);
D2(N,1) = -1;
D2(1,N) = -1;
D2 = D2 / (dx*dx);
% Construct P
P = zeros(N);
for i = 1:N
    if mod(i,2) == 1
        P((i + 1) / 2 , i) = 1;
    else
        P( N + 1 - i/2 ,i) = 1;
    end
end
PDP = P'*D2*P;
Pfx = P'*fx;

tic;
algo_mastery(PDP,Pfx); % Time algo mastery
tm = toc;

end