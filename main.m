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

% Initialize params
[robot, mpc, sim] = init_params();
ds = sim.ds;      % s-domain scale resolution
z  = zeros(mpc.z_dim, 1);      % initialize true error state
    
u  = zeros(3,1);     % initialize control input

% Generate reference path
ref = reference_trajectory(sim);

s = ref.s;      % s list, [0 : ds : S]
N = numel(s);   % loop N times

% Declare log
Z_log = zeros(mpc.z_dim, N);
U_log = zeros(mpc.u_dim, N);
t_log   = zeros(1, N);
X_log   = zeros(1, N);
Y_log   = zeros(1, N);
phi_log = zeros(1, N);

% Reference state, z_r, u_r
z_r = [0;
       0;
       ref.v_r;
       0;
       ref.phi_dot;
       0; 
       0]; 

u_r = zeros(3,1);

z(3) = ref.v_r;   

% linearize & discretize for prediction model in MPC
[Ad, Bd, Cd] = linearize_discretize(z_r, u_r, ds, ref);

t = 0;
t_last_u = 0;
control_period = mpc.Ts;   % = 0.1s


last_percent = 0;

% Start simulation
for s_idx = 1: N
    
    % ====== debug: nonlinear model input ======
    if any(~isfinite(z))
        fprintf('[ERROR] z corrupted before MPC at step %d\n', s_idx);
        disp(z.')
        break
    end

    % time increment in this step
    v_s   = z(3) * cos(z(2)) - z(4) * sin(z(2));  % same as f_continuous
    rho   = ref.rho;
    den   = rho - z(1);
    s_dot = rho * v_s / den;
    dt_ds = 1 / s_dot;
    dt    = dt_ds * ds;

    t = t + dt;      % accumulate time
    t_log(s_idx) = t;

    if (t - t_last_u >= control_period)

        % MPC
        u = mpc_solver(z, z_r, Ad, Bd, Cd, mpc);

        t_last_u = t;
    
        if any(~isfinite(u))
            fprintf('[ERROR] u is NaN at step %d\n', s_idx);
        end
    end

    % Update nonlinear Model from control input
    z_next = nonlinear_step(z, u, ds, ref);

    % ====== debug: nonlinear model output ======
    if any(~isfinite(z_next))
        fprintf('[ERROR] nonlinear_step produced NaN at step %d\n', s_idx);
        disp('z input:');
        disp(z.')
        disp('u input:');
        disp(u.')
        disp('z output:');
        disp(z_next.')
        break
    end

    z = z_next;  % Update the state for the next iteration

    

    % record
    Z_log(:, s_idx) = z;   % 7×N
    U_log(:, s_idx) = u;   % 3×N

    % transfer to global coordinate
    ref_k.s     = ref.s(s_idx);
    ref_k.rho   = ref.rho;      % scalar
    ref_k.phi_r = ref.phi_r(s_idx);
    [p, phi] = errorStateToWorld(z, ref_k);

    X_log(s_idx)   = p(1);
    Y_log(s_idx)   = p(2);
    % phi_log(s_idx) = phi;

    % --- Progress ---
    percent = floor(s_idx / N * 100);
    if percent ~= last_percent
        fprintf('\b\b\b%2d%%', percent);  % update progress
        last_percent = percent;
    end
end

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
%  Plot 4: all state 
% =========================================================================
figure;
plot(t_log, Z_log', 'LineWidth', 1.0);
xlabel('t [s]');
ylabel('state value');
title('All states (debug view)');
legend('e_y','e_\phi','v_x','v_y','\omega','a_x','a_y');
grid on;
% =========================================================================
%  Plot 5: global frame
% =========================================================================

figure;
% reference circle：
theta_ref = ref.phi_r - pi/2;        % for each sample
x_ref = ref.rho * cos(theta_ref);
y_ref = ref.rho * sin(theta_ref);

plot(x_ref, y_ref, 'w--', 'LineWidth', 1.5); hold on;
plot(X_log, Y_log, 'b--', 'LineWidth', 1.8);
axis equal;
xlabel('X [m]');
ylabel('Y [m]');
legend('Reference circle', 'Robot trajectory');
title('Path tracking in world frame');
grid on;
