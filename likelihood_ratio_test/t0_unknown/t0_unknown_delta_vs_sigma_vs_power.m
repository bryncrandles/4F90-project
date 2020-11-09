% t_0 unknown: 

mu_epsilon = 0;
sigma = 0:0.01:1;%0.001:1;
mu_0 = 10;
delta = 0:0.1:1.5;
obs_sigma_delta = zeros(length(sigma), length(delta), 1000);
T = 200;
t_0 = 100;

% loop to generate data for each combination of sigma and delta 
for i = 1:length(sigma)
    for j = 1:length(delta)
        S = zeros(1, 1000);
        for k = 1:length(S)
            X = zeros(1, T); % empty vector to store observed data 
            epsilon = normrnd(mu_epsilon, sigma(i), T, 1); % generate random normal error 
            for t = 1:T
                X(t) = mu_0 + delta(j)*heaviside(t - t_0) + epsilon(t);
            end
            var_hat_0 = var(X, 1);
            % create a vector of t_0 values, calculate mu_0_hat, delta_hat, var_hat_1
            t_0_val = 1:200; % vector of t_0 values to test
            % initialize vectors to store values of mu_0_hat, delta_hat, var_hat
            % and L (likelihood) that correspond to each t_0 value
            mu_0_hat_val = zeros(1, length(t_0_val));
            delta_hat_val = zeros(1, length(t_0_val));
            var_hat_1_val = zeros(1, length(t_0_val));
            L_val = zeros(1, length(t_0_val));
            for m = 1:length(t_0_val)
                mu_0_hat = mean(X(1:(t_0_val(m) - 1)));
                mu_0_hat_val(m) = mu_0_hat;
                delta_hat = mean(X(t_0_val(m):T)) - mean(X(1:(t_0_val(m) - 1)));
                delta_hat_val(m) = delta_hat;
                var_hat_1 = (1/T)*((t_0_val(m) - 1)*var(X(1:(t_0_val(m) - 1)), 1) + (T - t_0_val(m) + 1)*var(X(t_0_val(m):T), 1));
                var_hat_1_val(m) = var_hat_1;
                L_val(m) = (1/(2*pi*var_hat_1))^(T/2)*exp((-T/2));
            end
            if isinf(L_val(:))
               t_0_hat = find(var_hat_1 == min(var_hat_1));
            else
                t_0_hat = find(L_val == max(L_val));
            end
            if length(t_0_hat) > 1
                t_0_hat = t_0_hat(1);
            end
            var_hat_1 = var_hat_1_val(t_0_hat);
            % test statistic
            lambda_x = (var_hat_1/var_hat_0)^(T/2);
            chi_statistic = -2*log(lambda_x);
            S(k) = chi_statistic;
        end
        obs_sigma_delta(i, j, :) = S;
    end
end

%% Plot power of test as function of sigma and delta 
critical_val = chi2inv(0.95,1); % 3.84 
false_neg_rate_sigma_delta = sum(obs_sigma_delta < critical_val, 3)/size(obs_sigma_delta, 3);
power = 1 - false_neg_rate_sigma_delta;

figure()
p = pcolor(delta, sigma, power);
p.FaceColor = 'interp';
p.EdgeColor = 'none';
colorbar
title('t_0 Unknown: Delta vs Sigma vs Power', 'FontSize', 15)
xlabel('Delta', 'FontSize', 16)
ylabel('Sigma', 'FontSize', 16)

% save and exit 
saveas(p, 't_0_Unknown_Delta_vs_Sigma_vs_Power', 'png')

exit
