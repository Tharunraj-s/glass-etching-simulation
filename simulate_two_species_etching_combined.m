clc; clear;

% --- Parameters ---
f1 = 0.6; f2 = 0.4;
A1 = 1;   A2 = 2;
N = 300;                            % Number of layers
t_vals = [25, 64, 104, 144, 186];   % Time snapshots
vp = 1 / (f1/A1 + f2/A2);           % Vacancy front speed

fprintf('Theoretical wavefront speed vp = %.6f\n', vp);

% --- Initial Conditions ---
num_eqs = 2 * (N + 1);
Y0 = zeros(num_eqs, 1);
Y0(1) = f1;  % ψ1_0
Y0(2) = f2;  % ψ2_0

% --- ODE Definition ---
odefun = @(t, Y) compute_derivatives(Y, N, A1, A2, f1, f2);

% --- Solve ODE ---
tspan = [0 max(t_vals)];
[t_sol, Y_sol] = ode45(odefun, tspan, Y0);

% --- Plot ψ1_0(t) and ψ2_0(t) 
figure;
plot(t_sol, Y_sol(:,1), 'r-', 'LineWidth', 2); hold on;           % Species 1
plot(t_sol, Y_sol(:,2), 'b-', 'LineWidth', 2);                   % Species 2
xlabel('Time t');
ylabel('\psi_0(t)');
title('Species 1 and Species 2 at Layer 0');
legend('\psi^1_0(t)', '\psi^2_0(t)', 'Location', 'northeast');
grid on;

% --- Plot Vacant Sites vs Depth at Time Snapshots ---
n = 0:N;
figure; hold on;

% Define colors
col_total = 'k';        % black for total ψ_n
col_psi1  = [0, 0.45, 0.74];  % blue for ψ^1_n (e.g., RGB)
col_psi2  = [0.85, 0.33, 0.1]; % red-orange for ψ^2_n

for k = 1:length(t_vals)
    t_k = t_vals(k);
    Y_k = interp1(t_sol, Y_sol, t_k);

    psi1 = Y_k(1:2:end);
    psi2 = Y_k(2:2:end);
    psi  = psi1 + psi2;

    plot(n, psi, '-',  'Color', col_total, 'LineWidth', 1.5);     % ψ_n
    plot(n, psi1, '--', 'Color', col_psi1,  'LineWidth', 1.2);     % ψ^1_n
    plot(n, psi2, '-.', 'Color', col_psi2,  'LineWidth', 1.2);     % ψ^2_n
    xline(round(t_k * vp), 'k:');
end

xlabel('n (Depth into crystal)');
ylabel('\psi_n(t)');
title('Figure 3: Vacant Sites vs Depth');
legend('\psi_n(t)', '\psi^1_n(t)', '\psi^2_n(t)', 'Location', 'northeast');
xlim([0 N]);
ylim([0 0.08]);
grid on;


% --- Derivative function inline ---
function dYdt = compute_derivatives(Y, N, A1, A2, f1, f2)
    dYdt = zeros(size(Y));
    dYdt(1) = -A1 * Y(1);
    dYdt(2) = -A2 * Y(2);
    for n = 1:N
        i1 = 2*n + 1; i2 = 2*n + 2;
        ip1 = 2*(n-1) + 1; ip2 = 2*(n-1) + 2;
        flux = A1 * Y(ip1) + A2 * Y(ip2);
        dYdt(i1) = -A1 * Y(i1) + f1 * flux;
        dYdt(i2) = -A2 * Y(i2) + f2 * flux;
    end
end
