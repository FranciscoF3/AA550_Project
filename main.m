clear; clc; close all;

% =========================================================================
% Initialize params
% =========================================================================
[robot, mpc, sim] = init_params();
ds = sim.ds;      % s-domain scale resolution

% Generate reference path
ref = reference_trajectory(sim);

s = ref.s;      % s list, [0 : ds : S]
N = numel(s);   % loop N times

% Initialize robot state on reference
robot.x_cur  = ref.x_r(1);
robot.y_cur  = ref.y_r(1);
robot.phi_cur = ref.phi_r(1);

% Initialize true error state & control input
z  = zeros(mpc.z_dim, 1);      
z(3) = ref.v_r;             % v_x(0) cannot be zero!
u  = zeros(3,1);     

% Declare log
Z_log = zeros(mpc.z_dim, N);
U_log = zeros(mpc.u_dim, N);
t_log = zeros(1, N);
X_log = zeros(1, N);
Y_log = zeros(1, N);
D_log = zeros(3, N);    % disturbance log (3×N), like your groupmate

% Reference state z_r, u_r (for linearization and velocity reference)
z_r = [0;
       0;
       ref.v_r;
       0;
       ref.phi_dot;
       0; 
       0]; 

u_r = zeros(3,1);

% Linearize & discretize once (constant prediction model)
[Ad, Bd, Cd] = linearize_discretize(z_r, u_r, ds, ref);

% Terminal-state LQR weight
Q_lqr = diag([50 20 5 10 2 0.1 0.1]);
R_lqr = diag([0.01 0.01 0.10]);

[K_lqr, P_lqr, ~] = dlqr(Ad, Bd, Q_lqr, R_lqr);

%Using LQR to calculate W6

mpc.W6_ey   = P_lqr(1,1);
mpc.W6_ephi = P_lqr(2,2);

t = 0;
t_last_u = 0;
control_period = mpc.Ts;   % = 0.1s

last_percent = 0;

% ---------------- Disturbance initialization ----------------
d     = zeros(3,1);   % 3×1 disturbance: e.g. [d_vx; d_vy; d_w]
cfg_d = [];           % configuration struct if needed inside disturbance_step


% ---------- Pre-generate a shared disturbance profile ----------
rng(0);                      % fixed seed for reproducibility
dt_nom     = ds / sim.v_ref; % nominal time step for disturbance generation
d_gen      = zeros(3,1);
D_profile  = zeros(3, N);

for k = 1:N
    d_gen = disturbance_step(d_gen, dt_nom, cfg_d);
    D_profile(:,k) = d_gen;
end


% =========================================================================
% Start simulation
% =========================================================================
for s_idx = 1:N

    % time increment in this step (match f_continuous logic)
    v_s   = z(3) * cos(z(2)) - z(4) * sin(z(2));
    rho   = ref.rho;
    den   = rho - z(1);

    % simple protection against division by zero
    if abs(den) < 1e-3
        den = sign(den + 1e-3)*1e-3;
    end
    s_dot = rho * v_s / den;
    if abs(s_dot) < 1e-4
        s_dot = sign(s_dot + 1e-4)*1e-4;
    end

    dt_ds = 1 / s_dot;
    dt    = dt_ds * ds;

    t = t + dt;      % accumulate time
    t_log(s_idx) = t;

    % ------ MPC update at control_period ------
    if (t - t_last_u >= control_period)
        % store previous input for Δu term
        mpc.u_prev = u;

        % Quadprog-based MPC (unchanged signature)
        u = New_MPC_solver_QP(z, z_r, Ad, Bd, Cd, mpc);

        t_last_u = t;
    end

    % ------ Disturbance from pre-generated profile ------
    d = D_profile(:, s_idx);
    D_log(:, s_idx) = d;

    % ------ Nonlinear propagation with disturbance ------
    z_next = nonlinear_step(z, u, ds, ref, d);


    % Update the state for the next iteration
    z = z_next;

    % record
    Z_log(:, s_idx) = z;   % 7×N
    U_log(:, s_idx) = u;   % 3×N

    % Update robot state in world frame
    robot = robot_state_update(robot, z, ref, s_idx);
    X_log(s_idx) = robot.x_cur;
    Y_log(s_idx) = robot.y_cur;

    % --- Progress printout ---
    percent = floor(s_idx / N * 100);
    if percent ~= last_percent
        fprintf('\b\b\b%2d%%', percent);
        last_percent = percent;
    end
end

% =========================================================================
% Visualization
% =========================================================================
%  Plot 1: lateral & heading error
% =========================================================================
figure;
subplot(2,1,1);
plot(t_log, Z_log(1,:), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('e_y [m]');
title('Lateral error e_y(t)');
grid on;

subplot(2,1,2);
plot(t_log, Z_log(2,:), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('e_\phi [rad]');
title('Heading error e_\phi(t)');
grid on;

% =========================================================================
%  Plot 2: velocity states v_x, v_y, w
% =========================================================================
figure;
subplot(3,1,1);
plot(t_log, Z_log(3,:), 'LineWidth', 1.5); hold on;
yline(sim.v_ref, '--');
xlabel('Time [s]');
ylabel('v_x [m/s]');
title('Longitudinal velocity v_x(t)');
legend('v_x','v_x^{ref}');
grid on;

subplot(3,1,2);
plot(t_log, Z_log(4,:), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('v_y [m/s]');
title('Lateral velocity v_y(t)');
grid on;

subplot(3,1,3);
plot(t_log, Z_log(5,:), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('\omega [rad/s]');
title('Yaw rate \omega(t)');
grid on;

% =========================================================================
%  Plot 3: control inputs
% =========================================================================
figure;
subplot(3,1,1);
plot(t_log, U_log(1,:), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('\alpha [rad/s^3]');
title('Yaw jerk \alpha(t)');
grid on;

subplot(3,1,2);
plot(t_log, U_log(2,:), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('j_x [m/s^3]');
title('Longitudinal jerk j_x(t)');
grid on;

subplot(3,1,3);
plot(t_log, U_log(3,:), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('j_y [m/s^3]');
title('Lateral jerk j_y(t)');
grid on;

% =========================================================================
%  Plot 4: all states 
% =========================================================================
figure;
plot(t_log, Z_log', 'LineWidth', 1.0);
xlabel('t [s]');
ylabel('state value');
title('All states (debug view)');
legend('e_y','e_\phi','v_x','v_y','\omega','a_x','a_y');
grid on;

% =========================================================================
%  Plot 5: global frame path
% =========================================================================
figure;
theta_ref = ref.phi_r - pi/2;
x_ref = ref.rho * cos(theta_ref);
y_ref = ref.rho * sin(theta_ref);

plot(x_ref, y_ref, 'k--', 'LineWidth', 1.5); hold on;
plot(X_log, Y_log, 'b-', 'LineWidth', 1.8);

% start point
plot(X_log(1), Y_log(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
text(X_log(1), Y_log(1), '  Start', 'Color', 'g', 'FontSize', 12);

% end point
plot(X_log(end), Y_log(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(X_log(end), Y_log(end), '  End', 'Color', 'r', 'FontSize', 12);

axis equal;
xlabel('X [m]');
ylabel('Y [m]');
legend('Reference circle', 'Robot trajectory');
title('Path tracking in world frame');
grid on;

% =========================================================================
%  Plot 6: Disturbance vs time
% =========================================================================
figure;
plot(t_log, D_log', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('d');
legend('d_1','d_2','d_3');
title('Disturbance signals vs time');
grid on;

% ---- Save logs for comparison ----
mpc_results.t  = t_log;
mpc_results.Z  = Z_log;
mpc_results.X  = X_log;
mpc_results.Y  = Y_log;

save('mpc_results.mat','mpc_results');
