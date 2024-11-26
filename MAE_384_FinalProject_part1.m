% given parameters
h = 1; % time step 
T = 100; % total simulation time 
N = 1000; % total population
S0 = 990; % initial number of susceptible individuals
I0 = 10; % initial number of infected individuals
R0 = 0; % initial number of recovered individuals

% given parameters
parameters = {
    'Seasonal Influenza', 0.3, 0.1; % beta, gamma
    'COVID-19', 1, 0.1;            % beta, gamma
    'Measles', 2, 0.2              % beta, gamma
};

% derive sir model
sir_model = @(S, I, R, beta, gamma) [-beta*S*I/N; beta*S*I/N - gamma*I; gamma*I];

figure;
set(gcf, 'Position', [100, 100, 1200, 800]);

% parameter set loop
for p = 1:size(parameters, 1)
    disease = parameters{p, 1};
    beta = parameters{p, 2};
    gamma = parameters{p, 3};
    
    % initialize the given variables
    steps = T / h;
    S = zeros(steps+1, 1);
    I = zeros(steps+1, 1);
    R = zeros(steps+1, 1);
    t = (0:steps) * h;

    % set given initial conditions
    S(1) = S0;
    I(1) = I0;
    R(1) = R0;

    % Runge-Kutta 4th-order method
    for i = 1:steps
        k1 = h * sir_model(S(i), I(i), R(i), beta, gamma);
        k2 = h * sir_model(S(i) + 0.5*k1(1), I(i) + 0.5*k1(2), R(i) + 0.5*k1(3), beta, gamma);
        k3 = h * sir_model(S(i) + 0.5*k2(1), I(i) + 0.5*k2(2), R(i) + 0.5*k2(3), beta, gamma);
        k4 = h * sir_model(S(i) + k3(1), I(i) + k3(2), R(i) + k3(3), beta, gamma);

        S(i+1) = S(i) + (k1(1) + 2*k2(1) + 2*k3(1) + k4(1)) / 6;
        I(i+1) = I(i) + (k1(2) + 2*k2(2) + 2*k3(2) + k4(2)) / 6;
        R(i+1) = R(i) + (k1(3) + 2*k2(3) + 2*k3(3) + k4(3)) / 6;
    end

    % plot
    subplot(3, 1, p);
    plot(t, S, 'b', 'DisplayName', 'Susceptible (S)');
    hold on;
    plot(t, I, 'r', 'DisplayName', 'Infected (I)');
    plot(t, R, 'g', 'DisplayName', 'Recovered (R)');
    hold off;
    title(disease);
    xlabel('Time (days)');
    ylabel('Population');
    legend('Location', 'best');
    grid on;
end

sgtitle('SIR Model Simulation for Different Diseases');
