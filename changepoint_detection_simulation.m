% Simulate observed data X(t) = mu(t) + epsilon(t) where mu(t) = mu_0 +
% delta*H(t - t_0)

% Test hypothesis H0: delta = 0 vs H1: delta =/ 0 using likelihood ratio
% test
% set parameters
T = 200; % length of time (time -> infinity)
t_0 = 100; % time of change point corresponding to time
N = 1000; % number of simulations (N -> infinity)
mu_0 = 10; % mean 
mu_epsilon = 0; % mean of error
sigma = 1; % variance of error

% Test statistic of likelihood ratio test asympotically follows a Chi
% square distribution with 1 degree of freedom 
critical_val = chi2inv(0.95,1); % 3.84 

%% Repeat simulation with delta = 0 1000 times
delta = 0; % set delta = 0
obs_delta_zero = zeros(1, length(N)); % empty vector to store test statistics 

% loop to generate data for each length of time with delta = 0 and sigma = 1
for i = 1:N
    X = zeros(1, T); % empty vector to store observed data 
    epsilon = normrnd(mu_epsilon, sigma, 1, T); % generate random normal error 
    for t = 1:T
        X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
    end
    %var_0_hat = var(X, 1);
    %var_1_hat = (1/T)*((t_0 - 1)*var(X(1:t_0 - 1), 1) + (T - t_0 + 1)*var(X(t_0:T)));
    mu_0_hat = mean(X);
    xbar_1_hat = mean(X(1:t_0 - 1));
    xbar_2_hat = mean(X(t_0:T));
    S_0 = zeros(1, length(X));
    S_1 = zeros(1, (t_0 - 1));
    S_2 = zeros(1, (T - t_0 + 1));
    for m = 1:length(X)
        S_0(m) = (X(m) - mu_0_hat)^2;
    end
    for m = 1:(t_0 - 1)
        S_1(m) = (X(m) - xbar_1_hat)^2;
    end
    for m = t_0:T
        S_2(m - (t_0 - 1)) = (X(m) - xbar_2_hat)^2;
    end
    var_0_hat = sum(S_0);
    var_1_hat = sum(S_1) + sum(S_2);
    % test statistic
    lambda_x = (var_1_hat/var_0_hat)^(T/2);
    % -2*ln(lamda) approximately follows a Chi Square distribution with 1 df 
    chi_statistic = -2*log(lambda_x);
    obs_delta_zero(i) = chi_statistic;
end

%% Check false positive rate

false_pos_rate = sum(obs_delta_zero > critical_val)/length(obs_delta_zero);

% false positive rate = 0.03 with T = 200 and N = 1000

%% Repeat simulation 1000 times with delta = 0.01 to 10 using sigma = 1
T = 200;
t_0 = 100;
delta = 0.01:0.01:5;  % set delta vector
obs_delta = zeros(length(delta), 1000); % empty matrix to store data

%% loop to generate data with delta ranging from 0.01 to 10
for i = 1:length(delta)
    S = zeros(1, 1000); % empty vector to store chi statistics 
    for n = 1:1000 
        X = zeros(1, T); % empty vector to store observed data 
        epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error 
        for t = 1:T
            X(t) = mu_0 + delta(i)*heaviside(t - t_0) + epsilon(t);
        end
        %var_hat_0 = var(X, 1);
        %var_hat_1 = (1/T)*((t_0 - 1)*var(X(1:t_0 - 1), 1) + (T - t_0 + 1)*var(X(t_0:T), 1));
        mu_0_hat = mean(X);
        xbar_1_hat = mean(X(1:t_0 - 1));
        xbar_2_hat = mean(X(t_0:T));
        S_0 = zeros(1, length(X));
        S_1 = zeros(1, (t_0 - 1));
        S_2 = zeros(1, (T - t_0 + 1));
        for m = 1:length(X)
            S_0(m) = (X(m) - mu_0_hat)^2;
        end
        for m = 1:(t_0 - 1)
            S_1(m) = (X(m) - xbar_1_hat)^2;
        end
        for m = t_0:T
            S_2(m - (t_0 - 1)) = (X(m) - xbar_2_hat)^2;
        end
        var_0_hat = sum(S_0);
        var_1_hat = sum(S_1) + sum(S_2);
        % test statistic
        lambda_x = (var_1_hat/var_0_hat)^(T/2);
        chi_statistic = -2*log(lambda_x);
        S(n) = chi_statistic;
    end
    obs_delta(i, :) = S;
end

%% Plot power of the hypothesis test as function of delta
% (Recall power of a test: probability that test rejects H0 when H1 is true)
% power = 1 - false negative rate 

false_neg_rate_delta = sum(obs_delta < critical_val, 2)/length(obs_delta); 
power = 1 - false_neg_rate_delta;

figure()
plot(delta(1:100), power(1:100), 'b', 'LineWidth', 1.5) % plot from delta = 0 to 1 
title('Delta vs Power', 'FontSize', 14)
xlabel('Delta', 'FontSize', 16)
ylabel('Power', 'FontSize', 16)
set(gca,'FontSize', 15) % change font size of axis numbers

%% Power of the hypothesis test as a function of delta and sigma 
sigma = 0:0.001:1;
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
            mu_0_hat = mean(X);
            xbar_1_hat = mean(X(1:t_0 - 1));
            xbar_2_hat = mean(X(t_0:T));
            S_0 = zeros(1, length(X));
            S_1 = zeros(1, (t_0 - 1));
            S_2 = zeros(1, (T - t_0 + 1));
            for m = 1:length(X)
                S_0(m) = (X(m) - mu_0_hat)^2;
            end
            for m = 1:(t_0 - 1)
                S_1(m) = (X(m) - xbar_1_hat)^2;
            end
            for m = t_0:T
                S_2(m - (t_0 - 1)) = (X(m) - xbar_2_hat)^2;
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





    