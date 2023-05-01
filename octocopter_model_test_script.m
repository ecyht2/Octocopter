% Sourcing Drone Data
drone_data

% Setting values
X0 = zeros(12, 1);
X0(3, 1) = 4;
wh = sqrt(m * 9.81 / 8 / k_T);
U = wh * ones(8, 1);

% Running Simulation
sim("octocopter_model_test.slx");