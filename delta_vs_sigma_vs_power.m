%% Power of the hypothesis test as a function of delta and sigma 
mu_epsilon = 0;
sigma = 0:0.001:1;
mu_0 = 10;
delta = 0.01:0.1:1.5;
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
            %var_hat_0 = var(X, 1);
            %var_hat_1 = (1/T)*((t_0 - 1)*var(X(1:t_0 - 1), 1) + (T - t_0 + 1)*var(X(t_0:T), 1));
            x_0_hat = mean(X);
            x_1_hat = mean(X(1:t_0 - 1));
            x_2_hat = mean(X(t_0:T));
            S_0 = zeros(1, length(X));
            S_1 = zeros(1, (t_0 - 1));
            S_2 = zeros(1, (T - t_0 + 1));
            for m = 1:length(X)
                S_0(m) = (X(m) - x_0_hat)^2;
            end
            for m = 1:(t_0 - 1)
                S_1(m) = (X(m) - x_1_hat)^2;
            end
            for m = t_0:T
                S_2(m - (t_0 - 1)) = (X(m) - x_2_hat)^2;
            end
                var_0_hat = sum(S_0);
                var_1_hat = sum(S_1) + sum(S_2);
            % test statistic
            lambda_x = (var_hat_1/var_hat_0)^(T/2);
            chi_statistic = -2*log(lambda_x);
            S(k) = chi_statistic;
        end
            obs_sigma_delta(i, m, :) = S;
    end
end

%% Plot power of test as function of sigma and delta 
false_neg_rate_sigma_delta = sum(obs_sigma_delta < critical_val, 3)/size(obs_sigma_delta, 3);
power = 1 - false_neg_rate_sigma_delta;
%%
figure()
p = pcolor(delta, sigma, power);
p.FaceColor = 'interp';
p.EdgeColor = 'none';
colorbar
caxis = ('auto');
title('Delta vs Sigma vs Power', 'FontSize', 15)
xlabel('Delta', 'FontSize', 16)
ylabel('Sigma', 'FontSize', 16)





    