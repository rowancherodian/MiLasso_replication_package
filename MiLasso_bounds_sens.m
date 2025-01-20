%--------------------------------------------------------------------------
% MiLasso_bounds_sens
% Replication of the sensitivity analysis (supplement-Table 1) in 
% "Moran's I Lasso for models with spatially correlated data"
%--------------------------------------------------------------------------

clear variables
close all
clc

%--------------------------------------------------------------------------
% Set parameter values for sensitivity analysis
Z = 0.2;               % Fixed,       
alpha = 0.001;         % Required for o(.) behaviour limit as n -> infty

b_list = [1.5, 2, 5, 10];   % Sensitivity for 'b' parameter
n_list = [100, 200, 500];   % Sensitivity for sample size
sig_list = [0.5,1, 2, 4];   % Sensitivity for noise value

%--------------------------------------------------------------------------
% Preallocate output matrix and iterate over sensitivity dimensions
b_len = length(b_list);
n_len = length(n_list);
sig_len = length(sig_list);

bounds_mat = cell([sig_len*n_len+3,b_len+5]);
bounds_mat(:,:) = {repmat(' ',1,13)};

i0 = 1;
for k = 1:n_len                 % Iterate over values of n
    n = n_list(k);
    c_n = n^(-alpha);
    for j = 1:b_len             % Iterate over values of b
        b = b_list(j);
        for i = 1:sig_len       % Iterate over values of sigma
            sig = sig_list(i);

            bounds_mat(i0+i,1) = {sprintf('$\\sigma_v^2=%3g $ & ',sig)};

            L_bound = max(0,...
                log(2*b*sqrt(4*sig*log(n)/n) )/log(1/Z));   % Lower bound

            U_bound = c_n*log(n^(3/4))/log(1/Z);            % Upper bound

            bounds_mat(i0+i,2*j) = ...
                {sprintf('[%.3f,%.3f]',L_bound,U_bound)};

        end
        bounds_mat(:,2*j+1) = {' & '};

    end
    bounds_mat(i0,2) = {sprintf('&   $ n=%d $  ',n)};
    i0 = i0 + sig_len + 1;
end

%--------------------------------------------------------------------------
% Print table
fprintf('\\begin{tabular}{lcccc} \n')
fprintf('\\hline \n')

% Print header
for i = 1:b_len
    fprintf(' &  $ b=%g $  ',b_list(i))
end
fprintf('\\\\ \n')
fprintf('\\hline \n')

% Print table body
for i = 1:sig_len*n_len+3
    for j=1:b_len+4
        fprintf('%s',bounds_mat{i,j})
    end
    fprintf('\\\\ \n')
end
fprintf('\\hline \n')
fprintf('\\end{tabular} \n')
