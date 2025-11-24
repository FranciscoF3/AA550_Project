function ref = reference_trajectory(sim)
%   Generate a reference circle path for MPC with given s-resolution and laps.
% 
%   Args:
%       sim: (struct).
%           simulation config:
%               sim.ds    - (float), resolution in s-domain ([meter] per step)
%               sim.nLap  - (int)  , number of laps around the circle
%               sim.v_ref - (float), reference speed along the path (m/s)
% 
%   Returns:
%       ref: (struct).
%           A discretized reference trajectory. The path is defined as a circle 
%           with radius R in the global Cartesian frame. The trajectory is 
%           parameterized by the arc-length coordinate s, uniformly sampled 
%           with N points.
%               ref.s        - (1×N vector), arc-length samples along the trajectory.
%               ref.rho      - (float), circle path radius.
%               ref.phi_r    - (1×N vector), reference heading (tangent direction).
%               ref.v_r      - (float), reference forward speed profile.
%               ref.phi_dot  - (float), reference angular velocity, d(phi_r)/dt. 
    
    % Initialize params from sim config
    ds    = sim.ds;        % s-domain resolution (meter)
    n_lap = sim.n_lap;     % number of laps
    v_ref = sim.v_ref;     % reference speed along s
    R     = sim.rho;       % radius [meter]

    % total path length (multiple laps on the path)
    L_total = 2 * pi * R * n_lap;  

    % Discretize path
    s = 0 : ds : L_total;

    % Coordinate transformation: s-plane to global Cartesian coordinates
    theta = s / R;

    % Reference heading angle
    phi_r = theta + pi/2;
    phi_r = atan2(sin(phi_r), cos(phi_r));  % wrap to [-pi, pi]
    
    % Reference angular velocity
    phi_dot = v_ref / R;

    % Save into structure
    ref.s     = s;
    ref.rho   = R;
    ref.phi_r = phi_r;
    ref.v_r   = v_ref;
    ref.phi_dot = phi_dot;
end