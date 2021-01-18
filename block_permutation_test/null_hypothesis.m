% simulate AR(3) process

T = 200;
t_0 = 100; % time of change point 
delta = 0;
N = 1000; % number of simulations 
phi0 = 0;%5.5; % set to 5.5 for mean of 10
phi1 = 0.2; 
phi2 = 0.15;
phi3 = 0.1;
mu_epsilon = 0; % mean of error
sigma = 1; % variance of error
P = 10000;
stat = zeros(1, N);
alpha = 0.05;

for i = 1:N
    X = zeros(1, T);
    X(1) = 0;%10; 
    X(2) = 0;%10;
    X(3) = 0;%10;
    epsilon = normrnd(mu_epsilon, sigma, T, 1); % generate random normal error
    for t = 4:T
        X(t) = phi0 + delta*heaviside(t - t_0) + phi1*X(t - 1) + phi2*X(t - 2) + phi3*X(t - 3) + epsilon(t);
    end    
    xbar = mean(X);
    S = zeros(1, (T - 1));
    for t = 1:(T - 1)
        S(t) = sqrt(T/(t*(T - t)))*abs(sum(X(1:t) - xbar));
    end
    c_stat = max(S);
    c_stat_rand = zeros(1, P);
    c_stat_rand(1) = c_stat;
    for j = 2:P
        X_rand = X(randperm(T));
        S_N_rand = zeros(1, (T - 1));
        % loop to calculate S(t)
        for t = 1:(T - 1)
            S_N_rand(t) = sqrt(T/(t*(T - t)))*abs(sum(X_rand(1:t) - xbar));
        end
        c_stat_rand(j) = max(S_N_rand);
    end
    p_val = sum(c_stat_rand >= c_stat)/P;
    stat(i) = p_val;
end 

false_pos_rate = sum(stat < alpha)/N;
save('correlated_perm_test_valse_pos_rate.mat', 'false_pos_rate', '-mat')

exit 


    
