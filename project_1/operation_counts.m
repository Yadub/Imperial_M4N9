clear;

syms N k i;
Q_b = 4*(N + 1 - k);
S1 = symsum( Q_b , k, 1, N);
S1_lim = 4*N^2 / 2;
simplify(S1)

% test_case = 128;
% exact = symsum(S1,N,test_case,test_case)
% approx = symsum(S1_lim,N,test_case,test_case)
% (approx - exact)
% double ( (approx - exact) / exact )

R_b = 2*(N - i) + 1;
S2 = symsum( R_b , i, 1, N);
simplify(S2)
