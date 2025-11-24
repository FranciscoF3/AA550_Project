function [robot, mpc, sim] = init_params()
%   Initialize all robot, MPC, and simulation parameters used in the
%   robust MPC path-tracking framework for an omnidirectional mobile robot.
%
%   STRUCTURE:
%       robot : Physical and geometric parameters of the mecanum-wheel robot.
%       mpc   : MPC problem setup (state/input dimensions, horizons,
%               weighting matrices, constraints).
%       sim   : Simulation settings (reference speed, path resolution, laps).
%
%   Returns:
%       robot.mass     - Robot mass [kg]
%       robot.Jz       - Yaw inertia [kg·m^2]
%       robot.r        - Wheel radius [m]
%       robot.L, robot.l
%                     - Half-length/half-width of chassis [m]
%       robot.Fc, robot.Fv
%                     - Placeholder friction parameters
%
%       mpc.z_dim      - State dimension (7 for space-based model)
%       mpc.u_dim      - Input dimension (jerks + yaw accel)
%       mpc.sigma_dim  - Slack variable dimension for moving corridor
%       mpc.Ts         - MPC sampling period [s]
%       mpc.Kh         - Prediction horizon
%       mpc.Ch         - Control horizon
%       mpc.W1~W6      - Cost weights
%       mpc.u_min/max  - Input constraints
%       mpc.z_min/max  - State constraints
%
%       sim.ds         - Path resolution in s-domain
%       sim.v_ref      - Reference speed for space-based tracking
%       sim.n_lap      - Number of laps in closed-loop simulation
%
%   NOTE:
%       Modify this file when tuning control performance or adapting the
%       robot model. Keeping all parameters in one place prevents silent
%       inconsistencies across modules.

    % Robot parameters 
    robot.mass = 80;       % [kg]      
    robot.Jz   = 10;       % [kg m^2]  
    robot.r    = 0.1524;   % wheel radius [m]
    robot.L    = 0.241;    % half-length in x [m]
    robot.l    = 0.481;    % half-width in y [m]

    % Friction / motor model parameters (optional, TODO)
    robot.Fc = 0;          % Coulomb friction (placeholder)
    robot.Fv = 0;          % Viscous friction (placeholder)
    
    % ---------------------------------------------------------------------
    % MPC parameters
    % State & input 
    mpc.z_dim = 7;      % state dimenstion
    mpc.u_dim = 3;      % control input dimension
    
    % Moving corrider
    mpc.sigma_dim = 5; % moving corrider params

    % MPC sampling & horizon
    mpc.Ts = 0.1;   % control period[s]
    mpc.Kh = 15;    % prediction horizon 60
    mpc.Ch = 5;    % control horizon (<= Np) 10

    % Cost weights 
    mpc.W1 = 150;   % weight for lateral error ey
    mpc.W2 = 80;   % weight for heading error epsi
    mpc.W3 = diag([50 20 10]);
    mpc.W4 = eye(mpc.sigma_dim)*200;
    mpc.W5 = 1;
    mpc.W6 = 1;

    % Input constraints 
    mpc.u_min = [-1;       % alpha
                 -3;       % j_x
                 -3];      % j_y
    
    mpc.u_max = [ 1;       % alpha
                  3;       % j_x
                  3];      % j_y

    % State constraints
    mpc.z_min = [-0.5;      % e_y
                 -0.35;      % e_phi
                  0.1;      % v_x
                 -0.5;      % v_y
                 -1.5;      % w
                 -0.4;      % a_x
                 -0.4];     % a_y

    mpc.z_max = [ 0.5;      % e_y
                  0.35;      % e_phi
                  1.2;      % v_x
                  0.5;      % v_y
                  145;      % w??
                  0.4;      % a_x
                  0.4];     % a_y
    
    % ---------------------------------------------------------------------
    % Simulation parameters
    sim.ds    = 0.1;     % s-domain scale resolution (meter)
    sim.v_ref = 0.5;      % reference speed along s
    sim.n_lap = 1;        % number of laps
    sim.rho   = 5.0;      % radius [meter] 5 m
end