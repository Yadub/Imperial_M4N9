function [A, T] = optimized_house(A)

% This function performs the QR decomposition using the Householder reflections.  
% The input A is overwritten to give the matrix R.
% Input: A: complex mxn matrix with m>=n
% Outputs:  A: mxm upper-triangular matrix containing the coefficients to write
%             columns of A as a linear combination of the columns of Q (to be determined from W).
%           W: mxm lower triangular matrix containing the vectors v
%           defining the householder rreflections

[m, ~] = size(A);

e = zeros(m,1); % Offdiagonal elements
d = zeros(m,1); % To store the diagonal elements

for i = m:-1:2
    l = i-1;
    h = 0;
    scale = 0;
    if l > 1
        scale = sum( abs( A(i,1:l) ) );
        if scale == 0
            e(i) = A(i,i);
        else
            A(i,1:l) = A(i,1:l)/scale;
            h = sum( A(i,1:l).^2 );
            f = A(i,l);
            g = - sign(f) * sqrt(h);
            e(i) = scale * g;
            h = h - f * g;
            A(i,l) = f - g;
            f = 0;
            for j = 1:l
                A(j,l) = A(i,j) / h;
                g = sum( A(1:j,j)' .* A(i,1:j) );
                g = g + sum( A(j+1:l,j)' .* A(i,j+1:l) );
                e(j) = g / h;
                f = f + e(j)*A(i,j);
            end
            hh = f / (2.0*h);
            for j = 1:l
                f = A(i,j);
                g = e(j) - hh * f;
                e(j) = g;
                A(j,1:j) = A(j,1:j) - f * e(1:j)' - g * A(i,1:j);
            end
        end
    else
        e(i) = A(i,l);
        
    end
    d(i) = h;
end
d(1) = 0;
e(1) = 0;
for i = 1:m
    l = i-1;
    if d(i) ~= 0
        for j = 1:l
            g = sum( A(i,1:l)' .* A(1:l,j) );
            A(1:l,j) = A(1:l,j) - g * A(1:l,i);
        end
    end
    d(i) = A(i,i);
    A(i,i) = 1;
    A(i,1:l) = 0;
    A(1:l,i) = 0;
end

e2 = e(2:m);
T = diag(d) + diag(e2,1) + diag(e2,-1);

end 