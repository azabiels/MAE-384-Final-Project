% 384 Project Part 3 
%Initial Conditions 
beta = 0.3;
gamma = 0.1;
N = 1000;
S0 = 990;
I0 = 10;

% Time vector for 30 days
t = 1:30;

% non-linear SIR model to get true data I(t)
I_true = zeros(1, length(t));
I_true(1) = I0;

for i = 2:length(t)
    I_true(i) = I_true(i-1) * exp((beta*S0/N - gamma));
end

%Ln
ln_I_true = log(I_true);

% Linear least squares for 30 days
X = [ones(length(t), 1), t'];
Y = ln_I_true';
coefficients = (X' * X) \ (X' * Y);
ln_I0_30 = coefficients(1);
k_30 = coefficients(2);

% Estimate beta using k_30
beta_est_30 = (k_30 + gamma) * N / S0;

% Linear least squares for 10 days
t_10 = t(1:10);
ln_I_true_10 = ln_I_true(1:10);
X_10 = [ones(length(t_10), 1), t_10'];
coefficients_10 = (X_10' * X_10) \ (X_10' * Y(1:10));
ln_I0_10 = coefficients_10(1);
k_10 = coefficients_10(2);

% Estimate beta using k_10
beta_est_10 = (k_10 + gamma) * N / S0;

% Display results
fprintf('Estimated parameters using 30 days of data:\n');
fprintf('ln(I(0)): %.4f\n', ln_I0_30);
fprintf('k: %.4f\n', k_30);
fprintf('Estimated beta: %.4f\n', beta_est_30);

fprintf('\nEstimated parameters using 10 days of data:\n');
fprintf('ln(I(0)): %.4f\n', ln_I0_10);
fprintf('k: %.4f\n', k_10);
fprintf('Estimated beta: %.4f\n', beta_est_10);
