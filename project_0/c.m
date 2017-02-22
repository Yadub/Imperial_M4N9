% Yadu Bhageria
% CID: 00733164

m_vals = [128, 256, 512, 1024]; % Number of points
n = 9; % Number of polynmials to be computed

E_mmgs = zeros(3,4); % To store the error values for each polynomial with mmgs
E_smgs = zeros(3,4); % To store the error values for each polynomial with smgs

for index = 1:4
    m = m_vals(index);
    
    % Construct equally spaced points x_i
    delta_x = 2/m;
    x = zeros(m,1);
    for i = 1:m
        x(i) = delta_x / 2 + (i - 1) * delta_x - 1;
    end
    
    w = (1 - x.^2).^(-1/2); % Compute weight values for the error
    % Compute the exact values of the chebyshev polynomoials for current x
    T_m = zeros(m,n);
    T_m(:,1) = 1;
    T_m(:,2) = x;
    for i = 3:9
        T_m(:,i) = 2 * x .* T_m(:,i-1) - T_m(:,i-2);
    end
    
    Q_m = chebyshev_mmgs(x, n); % Compute the chebyshev polynomials using the modified MGS method
    
    % Construct equally spaced points over theta for x = cos(theta)
    for i = 1:m
        x(i) = cos( pi * ( 2 * i - 1) / ( 2 * m));
    end
    
    Q_s = chebyshev_smgs(x, n); % Compute the chebyshev polynomials using the standard MGS method
    
    % Compute the exact values of the chebyshev polynomoials for current x
    delta_theta = pi/m;
    T_s = zeros(m,n);
    T_s(:,1) = 1;
    T_s(:,2) = x;
    for i = 3:9
        T_s(:,i) = 2 * x .* T_s(:,i-1) - T_s(:,i-2);
    end
    
    poly_num = 2;
    for i = 1:3
        E_mmgs(i,index) = sqrt(delta_x * sum( ( T_m(:,poly_num + 1) - Q_m(:,poly_num + 1) ).^2 .* w));
        E_smgs(i,index) = sqrt(delta_theta *  sum(( T_s(:,poly_num + 1) - Q_s(:,poly_num + 1) ).^2));
        poly_num = poly_num*2;
    end
end

Table_mmgs = table(E_mmgs(:,1),E_mmgs(:,2),E_mmgs(:,3),E_mmgs(:,4), 'VariableNames',{'M128' 'M256' 'M512' 'M1024'}, 'RowNames',{'T_2' 'T_4' 'T_8'})
Table_smgs = table(E_smgs(:,1),E_smgs(:,2),E_smgs(:,3),E_smgs(:,4), 'VariableNames',{'M128' 'M256' 'M512' 'M1024'}, 'RowNames',{'T_2' 'T_4' 'T_8'})

% Plot results
clf;
figure;
for i = 1:3
    subplot(3,1,i);
    plot( m_vals, E_mmgs(1,:));
    title(['Errors using mmgs for T_' num2str(2^i)]);
end

figure;
for i = 1:3
    subplot(3,1,i);
    plot( m_vals, E_smgs(1,:));
    title(['Errors using smgs for T_' num2str(2^i)]);
end
