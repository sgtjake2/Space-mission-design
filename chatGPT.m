clc; clear; close all;

% Constants for the Earth-Moon system
mu = 0.012150585609624;  % Earth-Moon mass ratio
L2_x = 1.15535;          % Approximate L2 x-coordinate in normalized units

% Initial conditions using Richardson's third-order expansion
Az = 0.01;               % Chosen Z amplitude of the halo orbit
X0 = L2_x;               % Initial X position at L2
Y0 = 0;                  % Initial Y position
Z0 = Az;                 % Initial Z position
dX0 = 0;                 % Initial X velocity

% Compute initial Y velocity using analytical approximations
omega = sqrt(mu/(1-mu)); % Angular velocity
Vy0 = sqrt(omega^2 * (Z0^2));  % Approximate Y velocity
Vz0 = 0;                  % Initial Z velocity

% Initial state vector [X, Y, Z, dX, dY, dZ]
X_initial = [X0, Y0, Z0, dX0, Vy0, Vz0];

% Numerical integration of the orbit
tspan = [0, 4*pi];   % Time range
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t, state] = ode45(@(t, X) CR3BP_equations(t, X, mu), tspan, X_initial, options);

% Plot Halo Orbit
figure;
plot3(state(:,1), state(:,2), state(:,3), 'b', 'LineWidth', 1.5); hold on;
plot3(X0, Y0, Z0, 'ro', 'MarkerFaceColor', 'r'); % Initial point
xlabel('X (normalized)'); ylabel('Y (normalized)'); zlabel('Z (normalized)');
title('Halo Orbit around Earth-Moon L2');
grid on; axis equal;

% Function: CR3BP Equations of Motion
function dXdt = CR3BP_equations(~, X, mu)
    x = X(1); y = X(2); z = X(3);
    dx = X(4); dy = X(5); dz = X(6);
    
    r1 = sqrt((x + mu)^2 + y^2 + z^2); % Distance to Earth
    r2 = sqrt((x - (1-mu))^2 + y^2 + z^2); % Distance to Moon
    
    % Equations of motion
    ddx = 2*dy + x - (1-mu)*(x + mu)/r1^3 - mu*(x - (1-mu))/r2^3;
    ddy = -2*dx + y - (1-mu)*y/r1^3 - mu*y/r2^3;
    ddz = -(1-mu)*z/r1^3 - mu*z/r2^3;
    
    dXdt = [dx; dy; dz; ddx; ddy; ddz];
end
