clc; clear;

% --- Parameters ---
f1 = 0.3; f2 = 0.5; f3 = 0.2;
A1 = 1;   A2 = 2;   A3 = 3;
N = 300;                          % Number of layers
t_vals = [25, 64, 104, 144, 186]; % Time snapshots
vp = 1 / (f1/A1 + f2/A2 + f3/A3); % Wavefront speed

fprintf('Theoretical wavefront speed vp = %.6f\n', vp);

% --- Initial Conditions ---
num_eqs = 3 * (N + 1);  % 3 species × (N+1 layers)
Y0 = zeros(num_eqs, 1);
Y0(1) = f1;  % ψ1_0
Y0(2) = f2;  % ψ2_0
Y0(3) = f3;  % ψ3_0

% --- ODE Definition ---
odefun = @(t, Y) compute_derivatives_3species(Y, N, A1, A2, A3, f1, f2, f3);

% --- Solve ODE ---
tspan = [0 max(t_vals)];
[t_sol, Y_sol] = ode45(odefun, tspan, Y0);

% --- Plot ψ1_0(t), ψ2_0(t), ψ3_0(t) ---
figure;
plot(t_sol, Y_sol(:,1), 'r-', 'LineWidth', 2); hold on;
plot(t_sol, Y_sol(:,2), 'b-', 'LineWidth', 2);
plot(t_sol, Y_sol(:,3), 'g-', 'LineWidth', 2);
xlabel('Time t'); ylabel('\psi^i_0(t)');
title('Species Fractions at Surface Layer n = 0');
legend('\psi^1_0(t)', '\psi^2_0(t)', '\psi^3_0(t)');
grid on;

% --- Plot Vacant Sites vs Depth at Time Snapshots ---
n = 0:N;
figure; hold on;

for k = 1:length(t_vals)
    t_k = t_vals(k);
    Y_k = interp1(t_sol, Y_sol, t_k);

    psi1 = Y_k(1:3:end);
    psi2 = Y_k(2:3:end);
    psi3 = Y_k(3:3:end);
    psi  = psi1 + psi2 + psi3;

    plot(n, psi, 'k-', 'LineWidth', 1.5);        % ψ_n
    plot(n, psi1, 'r--', 'LineWidth', 1.2);      % ψ1_n
    plot(n, psi2, 'b-.', 'LineWidth', 1.2);      % ψ2_n
    plot(n, psi3, 'g:',  'LineWidth', 1.5);      % ψ3_n
    xline(round(t_k * vp), 'k:');                % wavefront
end

xlabel('n (Depth into crystal)');
ylabel('\psi_n(t)');
title('Vacant Sites vs Depth');
legend('\psi_n(t)', '\psi^1_n(t)', '\psi^2_n(t)', '\psi^3_n(t)', 'Location', 'northeast');
xlim([0 N]);
ylim([0 0.08]);
grid on;

% --- Derivative Function for 3 Species ---
function dYdt = compute_derivatives_3species(Y, N, A1, A2, A3, f1, f2, f3)
    dYdt = zeros(size(Y));

    % Layer 0
    dYdt(1) = -A1 * Y(1);  % ψ1_0
    dYdt(2) = -A2 * Y(2);  % ψ2_0
    dYdt(3) = -A3 * Y(3);  % ψ3_0

    % Layers 1 to N
    for n = 1:N
        i1 = 3*n + 1; i2 = 3*n + 2; i3 = 3*n + 3;
        ip1 = 3*(n-1) + 1; ip2 = 3*(n-1) + 2; ip3 = 3*(n-1) + 3;

        flux = A1 * Y(ip1) + A2 * Y(ip2) + A3 * Y(ip3);

        dYdt(i1) = -A1 * Y(i1) + f1 * flux;
        dYdt(i2) = -A2 * Y(i2) + f2 * flux;
        dYdt(i3) = -A3 * Y(i3) + f3 * flux;
    end
end
