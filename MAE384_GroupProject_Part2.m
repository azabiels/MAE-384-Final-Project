clear all; clc;

% Parameters
T = 100; % Total run time in days
h_coarse = 2; % Coarse time step
h_fine = 1; % Fine time step
beta = 0.3; % Infection rate
gamma = 0.1; % Recovery rate

% Initial Conditions
S0 = 990;
I0 = 10;
R0 = 0; % Initial populations
N = S0 + I0 + R0; % Total population

% SIR Model Differential Equations
dS_dt = @(S, I) -beta * S .* I / N;
dI_dt = @(S, I) beta * S .* I / N - gamma * I;
dR_dt = @(I) gamma * I;

% Runge-Kutta 4th Order Solver
function [S, I, R, t] = RK4_solver(S0, I0, R0, h, T, dS_dt, dI_dt, dR_dt)
    num_steps = T / h;
    t = 0:h:T; % Time array
    S = zeros(1, num_steps + 1); I = zeros(1, num_steps + 1); R = zeros(1, num_steps + 1);
    S(1) = S0; I(1) = I0; R(1) = R0;
    
    for n = 1:num_steps
        k1_S = h * dS_dt(S(n), I(n));
        k1_I = h * dI_dt(S(n), I(n));
        k1_R = h * dR_dt(I(n));
        
        k2_S = h * dS_dt(S(n) + 0.5 * k1_S, I(n) + 0.5 * k1_I);
        k2_I = h * dI_dt(S(n) + 0.5 * k1_S, I(n) + 0.5 * k1_I);
        k2_R = h * dR_dt(I(n) + 0.5 * k1_I);
        
        k3_S = h * dS_dt(S(n) + 0.5 * k2_S, I(n) + 0.5 * k2_I);
        k3_I = h * dI_dt(S(n) + 0.5 * k2_S, I(n) + 0.5 * k2_I);
        k3_R = h * dR_dt(I(n) + 0.5 * k2_I);
        
        k4_S = h * dS_dt(S(n) + k3_S, I(n) + k3_I);
        k4_I = h * dI_dt(S(n) + k3_S, I(n) + k3_I);
        k4_R = h * dR_dt(I(n) + k3_I);
        
        S(n+1) = S(n) + (k1_S + 2*k2_S + 2*k3_S + k4_S) / 6;
        I(n+1) = I(n) + (k1_I + 2*k2_I + 2*k3_I + k4_I) / 6;
        R(n+1) = R(n) + (k1_R + 2*k2_R + 2*k3_R + k4_R) / 6;
    end
end

% Baseline Case for h=1 and h=2
[S_fine, I_fine, R_fine, t_fine] = RK4_solver(S0, I0, R0, h_fine, T, dS_dt, dI_dt, dR_dt);
[S_coarse, I_coarse, R_coarse, t_coarse] = RK4_solver(S0, I0, R0, h_coarse, T, dS_dt, dI_dt, dR_dt);

% Newton's Interpolation
odd_days = 1:2:T; % Days to interpolate
f_linear = @(y0, y1, x0, x1, x) y0 + (y1 - y0) * (x - x0) / (x1 - x0);
f_quadratic = @(y0, y1, y2, x0, x1, x2, x) y0 + (x - x0) * ((y1 - y0) / (x1 - x0)) + (x - x0) * (x - x1) * (((y2 - y1) / (x2 - x1) - (y1 - y0) / (x1 - x0)) / (x2 - x0));

% Linear and Quadratic Interpolation
S_interp_lin = zeros(1, length(odd_days));
S_interp_quad = zeros(1, length(odd_days));
I_interp_lin = zeros(1, length(odd_days));
I_interp_quad = zeros(1, length(odd_days));
R_interp_lin = zeros(1, length(odd_days));
R_interp_quad = zeros(1, length(odd_days));

for k = 1:length(odd_days)
    idx = find(t_coarse <= odd_days(k), 1, 'last');
    
    % Linear Interpolation (safe for all cases)
    S_interp_lin(k) = f_linear(S_coarse(idx), S_coarse(idx+1), t_coarse(idx), t_coarse(idx+1), odd_days(k));
    I_interp_lin(k) = f_linear(I_coarse(idx), I_coarse(idx+1), t_coarse(idx), t_coarse(idx+1), odd_days(k));
    R_interp_lin(k) = f_linear(R_coarse(idx), R_coarse(idx+1), t_coarse(idx), t_coarse(idx+1), odd_days(k));
    
    % Quadratic Interpolation (only if idx+2 exists)
    if idx + 2 <= length(t_coarse)
        S_interp_quad(k) = f_quadratic(S_coarse(idx), S_coarse(idx+1), S_coarse(idx+2), t_coarse(idx), t_coarse(idx+1), t_coarse(idx+2), odd_days(k));
        I_interp_quad(k) = f_quadratic(I_coarse(idx), I_coarse(idx+1), I_coarse(idx+2), t_coarse(idx), t_coarse(idx+1), t_coarse(idx+2), odd_days(k));
        R_interp_quad(k) = f_quadratic(R_coarse(idx), R_coarse(idx+1), R_coarse(idx+2), t_coarse(idx), t_coarse(idx+1), t_coarse(idx+2), odd_days(k));
    else
        % Use linear interpolation as a fallback
        S_interp_quad(k) = S_interp_lin(k);
        I_interp_quad(k) = I_interp_lin(k);
        R_interp_quad(k) = R_interp_lin(k);
    end
end

% Compute L2 Errors
L2_error = @(V_interp, V_model) sqrt(sum((V_interp - V_model).^2) / length(V_interp));
S_error_lin = L2_error(S_interp_lin, S_fine(odd_days));
S_error_quad = L2_error(S_interp_quad, S_fine(odd_days));
I_error_lin = L2_error(I_interp_lin, I_fine(odd_days));
I_error_quad = L2_error(I_interp_quad, I_fine(odd_days));
R_error_lin = L2_error(R_interp_lin, R_fine(odd_days));
R_error_quad = L2_error(R_interp_quad, R_fine(odd_days));

error_table = table(["S_Error"; "I_Error"; "R_Error"], [S_error_lin; I_error_lin; R_error_lin], [S_error_quad; I_error_quad; R_error_quad], 'VariableNames', {'Population', 'Linear', 'Quadratic'});
disp(error_table);