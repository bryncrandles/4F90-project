%% Power of the hypothesis test as a function of delta and sigma 
T = 200;
t_0 = 100;
mu_epsilon = 0;
sigma = 0:0.01:1;%0.001:1;
mu_0 = 10;
delta = 0:0.1:1.5;%0.01:1.5;
obs_sigma_delta = zeros(length(sigma), length(delta), 1000);
N = 1000;

% Test statistic of likelihood ratio test asympotically follows a Chi
% square distribution with 1 degree of freedom 
critical_val = chi2inv(0.95,1); % 3.84 

% loop to generate data for each combination of sigma and delta 
for i = 1:length(sigma)
    for j = 1:length(delta)
        S = zeros(1, N);
        for k = 1:length(S)
            X = zeros(1, T); % empty vector to store observed data 
            epsilon = normrnd(mu_epsilon, sigma(i), T, 1); % generate random normal error 
            for t = 1:T
                X(t) = mu_0 + delta(j)*heaviside(t - t_0) + epsilon(t);
            end
            var_0_hat = var(X, 1);
            var_1_hat = (1/T)*((t_0 - 1)*var(X(1:(t_0 - 1)), 1) + (T - t_0 + 1)*var(X(t_0:T), 1));
            % test statistic
            lambda_x = (var_1_hat/var_0_hat)^(T/2);
            chi_statistic = -2*log(lambda_x);
            S(k) = chi_statistic;
        end
            obs_sigma_delta(i, j, :) = S;
    end
end

% Plot power of test as function of sigma and delta 
false_neg_rate_sigma_delta = sum(obs_sigma_delta < critical_val, 3)/size(obs_sigma_delta, 3);
power = 1 - false_neg_rate_sigma_delta;
figure()
p = pcolor(delta, sigma, power);
p.FaceColor = 'interp';
p.EdgeColor = 'none';
colorbar
title('Delta vs Sigma vs Power', 'FontSize', 15)
xlabel('Delta', 'FontSize', 16)
ylabel('Sigma', 'FontSize', 16)

% save figure
saveas(p, 'Delta_vs_Sigma_vs_Power', 'png')

exit


    