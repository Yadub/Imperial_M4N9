clear

% This part of the code sets up the points where the monomials in the
% Vandermonde matrix will be evaluate

M = 1024; %number of points (number of rows of A)

% These are the equispace points in x used with the weighted norm. 
dx = 2/M;
xw = (-1 + dx/2:dx:1 - dx/2)';

% These are the equispace points in theta used with the standard MGS
dtheta = pi/M;
theta = (pi-dtheta/2:-dtheta:dtheta/2)';
x = cos(theta); 

%------------------------------------------------------------%

% The next few lines set up the Vandermonde matrix 

N = 9; % The number of monomials (number of colums of A).  The highest degree polynomial is N-1.
A = zeros(M,N); % Allocate memory for A
Aw = zeros(M,N);


for j=1:N
    A(:,j) = x.^(j-1); 
    Aw(:,j) = xw.^(j-1); 
end

%------------------------------------------------------------%

% Perform QR using modified Gram-Schmidt
[Qm, Rm] = mgs(A);

% Sets up the weighting matrix for the discrete weighted norm.
W = 1./sqrt(1 - xw.^2);
% Perform QR using modified Gram-Schmidt with weighted norms.
[Qw, Rw] = mgs_wnorm(Aw,W);

% This part of the code uses the recursion relation to obtain the values of
% the Chebyshev polynomial of degree N-1 at both sets of points.

T0 = ones(M,1);
T1 = x;

T0w = ones(M,1);
T1w = xw;

for i=3:N
    T2 = 2.*x.*T1 - T0;
    T0 = T1;
    T1 = T2;
    
    T2w = 2.*xw.*T1w - T0w;
    T0w = T1w;
    T1w = T2w;
end

% Here we compute the errors making sure to appropriately normalize.

Tnum = T2(end)*Qm(:,N)/Qm(end,N);
error1 = sqrt(dtheta*sum((T2 - Tnum).^2))

Tnum = T2w(end)*Qw(:,N)/Qw(end,N);
error2 = sqrt(dx)*wnorm(T2w - Tnum,W)
