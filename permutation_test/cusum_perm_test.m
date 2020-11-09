%% Verify the false positive rate of permutation test for cusum statistic

% Generate time series data 
T = 200; % length of time 
t_0 = 100; % time of change point 
mu_0 = 10;
N = 1000; % number of simulations 
mu_epsilon = 0; % mean of error
sigma = 1; % variance of error
delta = 0;
P = 10000;
stat = zeros(1, N);
alpha = 0.05;

for i = 1:N
    X = zeros(1, T);
    epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error
    for t = 1:T
        X(t) = mu_0 + delta*heaviside(t - t_0) + epsilon(t);
    end
    xbar = mean(X);
    S = zeros(1, (T - 1));
    for t = 1:(T - 1)
        S(t) = sqrt(T/(t*(T - t)))*abs(sum(X(1:t) - xbar));
    end
    c_stat = max(S);
    c_stat_rand = zeros(1, P);
    for j = 1:P
        X_rand = X(randperm(T));
        S_N_rand = zeros(1, (T - 1));
        % loop to calculate S(t)
        for t = 1:(T - 1)
            S_N_rand(t) = sqrt(T/(t*(T - t)))*abs(sum(X_rand(1:t) - xbar));
        end
        c_stat_rand(j) = max(S_N_rand);
    end
    p_val = sum(c_stat_rand > c_stat)/P;
    stat(i) = p_val;
end 

false_pos_rate = sum(stat < alpha)/N;
save('perm_test_valse_pos_rate', 'false_pos_rate', '-mat')

% when N = 100, P = 10000, false positive rate is 0.04.

%% Power of the test as delta ranges from 0 to 5

delta = 0:0.1:5;
obs_delta = zeros(length(delta), N);

for i = 1:length(delta)
    stat = zeros(1, N);
    for j = 1:N
        X = zeros(1, T);
        epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error
        for t = 1:T
            X(t) = mu_0 + delta(i)*heaviside(t - t_0) + epsilon(t);
        end
        xbar = mean(X);
        S = zeros(1, (T - 1));
        for t = 1:(T - 1)
            S(t) = sqrt(T/(t*(T - t)))*abs(sum(X(1:t) - xbar));
        end
        c_stat = max(S);
        c_stat_rand = zeros(1, P);
        for k = 1:P
            X_rand = X(randperm(T));
            S_N_rand = zeros(1, (T - 1));
            % loop to calculate S(t)
            for t = 1:(T - 1)
                S_N_rand(t) = sqrt(T/(t*(T - t)))*abs(sum(X_rand(1:t) - xbar));
            end
            c_stat_rand(k) = max(S_N_rand);
        end
        p_val = sum(c_stat_rand > c_stat)/P;
        stat(j) = p_val;
    end 
    obs_delta(i, :) = stat;
end

%% plot power 
false_neg_rate = sum(obs_delta > alpha, 2)/N;
power = 1 - false_neg_rate;

plot = plot(delta(1:100), power(1:100), 'b', 'LineWidth', 1.5); % plot from delta = 0 to 1 
title('Permutation Test: Delta vs Power', 'FontSize', 14)
xlabel('Delta', 'FontSize', 16)
ylabel('Power', 'FontSize', 16)
set(gca,'FontSize', 15) % change font size of axis numbers
saveas(plot, 'perm_test_delta_vs_power', 'png')



