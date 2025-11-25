% Define state and input:
%   z:  7x1 matrix. 
%       Error-state vector:
%       [e_y; e_phi; v_x; v_y; w; a_x; a_y]^T
%           e_y   : lateral error
%           e_phi : heading error
%           v_x   : longitudinal velocity in robot frame
%           v_y   : lateral velocity in robot frame
%           w     : yaw rate
%           a_x   : longitudinal acceleration
%           a_y   : lateral acceleration
% 
%   u:  3x1 matrix. 
%       Control input at current step:
%       [alpha; j_x; j_y]^T
%           alpha : yaw jerk 
%           j_x   : longitudinal jerk
%           j_y   : lateral jerk
% 
% =========================================================================

clear; clc;
% Initialize params
[robot, mpc, sim] = init_params();
ds = sim.ds;      % s-domain scale resolution

% Generate reference path
ref = reference_trajectory(sim);

s = ref.s;      % s list, [0 : ds : S]
N = numel(s);   % loop N times

% Initialize robot state
robot.x_cur = ref.x_r(1);
robot.y_cur = ref.y_r(1);
robot.phi_cur = ref.phi_r(1);

% Initialize true error state & control input
z  = zeros(mpc.z_dim, 1);      
z(3) = ref.v_r;             % v_x(0) cannot be zero!
u  = zeros(3,1);     

% Initialize recorder
% vis = visualizer_init(ref, N);

% Declare log
Z_log = zeros(mpc.z_dim, N);
U_log = zeros(mpc.u_dim, N);
t_log = zeros(1, N);
X_log = zeros(1, N);
Y_log = zeros(1, N);
D_log = zeros(3, N);  

% Reference state, z_r, u_r
z_r = [0;
       0;
       ref.v_r;
       0;
       ref.phi_dot;
       0; 
       0]; 

u_r = zeros(3,1);

% linearize & discretize for prediction model in MPC
[Ad, Bd, Cd] = linearize_discretize(z_r, u_r, ds, ref);

t = 0;
t_last_u = 0;
control_period = mpc.Ts;   % = 0.1s

last_percent = 0;

% Disturbance
d = zeros(3,1);       % 初始 disturbance
cfg_d = [];           % 或 cfg_d = struct(...)，如果你要自訂參數

% Start simulation
for s_idx = 1: N
    
    % ====== debug: nonlinear model input ======
    % if any(~isfinite(z))
    %     fprintf('[ERROR] z corrupted before MPC at step %d\n', s_idx);
    %     disp(z.')
    %     break
    % end

    % time increment in this step
    v_s   = z(3) * cos(z(2)) - z(4) * sin(z(2));  
    rho   = ref.rho;
    den   = rho - z(1);
    s_dot = rho * v_s / den;

    eps_s = 1e-6;
    if abs(s_dot) < eps_s
        s_dot = sign(s_dot + eps_s) * eps_s;
    end

    dt_ds = 1 / s_dot;
    dt    = dt_ds * ds;

    t = t + dt;      % accumulate time
    t_log(s_idx) = t;

    if (t - t_last_u >= control_period)
        
        tic;
        % MPC
        u = mpc_solver(z, z_r, Ad, Bd, Cd, mpc, ref);
        t_solve = toc;

        t_last_u = t;
    
        % if any(~isfinite(u))
        %     fprintf('[ERROR] u is NaN at step %d\n', s_idx);
        % end
    end
    
    % disturbance at this step
    d = disturbance_step(d, dt, cfg_d);  % 3x1
    D_log(:, s_idx) = d;

    % Update nonlinear Model from control input
    % z_next = nonlinear_step(z, u, ds, ref);
    z_next = nonlinear_step(z, u, ds, ref, d);

    % % ====== debug: nonlinear model output ======
    % if any(~isfinite(z_next))
    %     fprintf('[ERROR] nonlinear_step produced NaN at step %d\n', s_idx);
    %     disp('z input:');
    %     disp(z.')
    %     disp('u input:');
    %     disp(u.')
    %     disp('z output:');
    %     disp(z_next.')
    %     break
    % end

    % record
    Z_log(:, s_idx) = z;   % 7×N
    U_log(:, s_idx) = u;   % 3×N

    % Update robot state
    robot = robot_state_update(robot, z, ref, s_idx);
    X_log(s_idx) = robot.x_cur;
    Y_log(s_idx) = robot.y_cur;

    % ==== Real-time visualization (packed into function) ====
    % [vis, stop_sim] = visualizer_update(vis, s_idx, t, robot, z, u);
    % 
    % if stop_sim
    %     fprintf('\n⚠️  Simulation aborted at step %d due to abnormal behavior.\n', s_idx);
    %     break;
    % end

    % --- Progress ---
    percent = floor(s_idx / N * 100);
    if percent ~= last_percent
        fprintf('\b\b\b%2d%%', percent);  % update progress
        last_percent = percent;
    end
    
    % Update the state for the next iteration
    z = z_next;  
end

% Visualization
%% ========================================================================
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

%% ========================================================================
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

%% ========================================================================
%  Plot 3: Control inputs v.s. time
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


%% ========================================================================
%  Plot 4: Error state v.s. time
% =========================================================================
figure;
plot(t_log, Z_log', 'LineWidth', 1.0);
xlabel('t [s]');
ylabel('state value');
title('All states (debug view)');
legend('e_y','e_\phi','v_x','v_y','\omega','a_x','a_y');
grid on;

%% =========================================================================
%  Plot 5: global frame
% =========================================================================

figure;
% reference circle：
theta_ref = ref.phi_r - pi/2;        % for each sample
x_ref = ref.rho * cos(theta_ref);
y_ref = ref.rho * sin(theta_ref);

plot(x_ref, y_ref, 'w--', 'LineWidth', 1.5); hold on;
plot(X_log, Y_log, 'b--', 'LineWidth', 1.8);

% start point
plot(X_log(1), Y_log(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
text(X_log(1), Y_log(1), '  Start', 'Color', 'g', 'FontSize', 12);

% end point
plot(X_log(end), Y_log(end), 'rx', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(X_log(end), Y_log(end), '  End', 'Color', 'r', 'FontSize', 12);

% % plot obstacle
% % 取 obstacle 參數
% cx = mpc.obs.cx;
% cy = mpc.obs.cy;
% R  = mpc.obs.R_safe;
% 
% % 畫 obstacle center（用 x 標記）
% plot(cx, cy, 'rx', 'LineWidth', 2, 'MarkerSize', 10);
% 
% % 畫安全範圍圓
% theta = linspace(0, 2*pi, 200);
% x_circle = cx + R * cos(theta);
% y_circle = cy + R * sin(theta);
% plot(x_circle, y_circle, 'r-', 'LineWidth', 1.5);
% 
% % optional: 填滿淡紅色區域（更顯眼）
% patch(x_circle, y_circle, 'r', ...
%       'FaceAlpha', 0.1, ...   % 透明度
%       'EdgeColor', 'none');

legend('Reference path', 'Robot path', 'Obstacle center', 'Safe region');
xlabel('x [m]');
ylabel('y [m]');
title('Path tracking with obstacle and safe region');

axis equal;
% xlabel('X [m]');
% ylabel('Y [m]');
% legend('Reference circle', 'Robot trajectory');
% title('Path tracking in world frame');
grid on;

%% ========================================================================
%  Plot 6: Disturbance v.s. time
% =========================================================================
figure;
plot(t_log, D_log', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('d_{ay}');
legend('lateral dis', 'longitude dis', 'angular dis');
title('Disturbance');
grid on;
