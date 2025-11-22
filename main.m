close all; clear; clc;

% Initialize params
[mpc, sim] = init_params();
ds = sim.ds;      % s-domain scale resolution
z  = sim.z0;      % initialize true error state

% Generate reference path
ref = reference_trajectory(sim);

s = ref.s;      % s list
N = numel(s);   % loop N times

% log
Z_log = zeros(sim.z_dim, N);
U_log = zeros(sim.u_dim, N);

% Start sim
for s_idx = 1: N
    % Sampling
    target_ref = reference_sampling(ref, s_idx);

    % z_r, u_r
    z_r = [0;
           0;
           ref.v_r;
           0;
           0;
           0; 
           ref.phi_dot; ];  % ??

    u_r = [0; 0; 0;];
    
    % linearize & discretize for prediction model in MPC
    [Ad, Bd, Cd] = linearize_discretize(z_r, u_r, ds, ref);

    % MPC

    % Update nonlinear Model
    for i = 1: mpc.Ch
        z = nonlin_model(u(i), z, delta_s, ref);

        % record
        Z_log(k,:) = z.';
        U_log(k,:) = u.';
    end
end

% Visualization
plot_results(s, Z_log, U_log, ref, sim, mpc, robot);
