function z_next = nonlinear_step(z, u, delta_s, ref, d)
%/Update error state z(s) using the nonlinear space-based model.
% 
% Args:
%   u_cur:  3x1 matrix. 
%           Control input at current step:
%           [alpha; j_x; j_y]^T
% 
%   z_cur:  7x1 matrix. 
%           Error-state vector in path-aligned frame:
%           [e_y; e_phi; v_x; v_y; w; a_x; a_y]^T
%          
%   delta_s: scalar.
%            Path step length Δs used for space-based integration.
% 
%   ref: struct.
%        Reference path information at current step, e.g.:
%           ref.phi_dot : reference heading rate (φ̇_r)
%           ref.rho     : path parameter (e.g. curvature-related term in ṡ definition)
% 
% Returns:
%   true_z_next: 7x1 matrix. 
%                Updated error-state at the next path step s + Δs:
%                [e_y; e_phi; v_x; a_x; v_y; a_y; w]^T
% /% 
    
    % Whether include d in inputs
    if nargin < 5
        d = [0; 0; 0];
    end

    % Nonlinear continuous model
    dz_ds = f_continuous(z, u, ref, d);  

    % Update z(s)
    z_next = z + dz_ds * delta_s;   % Euler integration
   
end



function dz_ds = f_continuous_d(z, u, ref, d)
%   Computes space-based nonlinear dynamics dz/ds = f(z,u) with disturbance
%   
%   Args:
%       u_cur:  3x1 matrix. 
%               Control input at current step:
%               [alpha; j_x; j_y]^T
% 
%       z_cur:  7x1 matrix. 
%               Error-state vector in path-aligned frame:
%               [e_y; e_phi; v_x; v_y; w; a_x; a_y]^T
%       ref:    struct.
%               Reference path information at current step, e.g.:
%                   ref.phi_dot : reference heading rate (φ̇_r)
%                   ref.rho     : path radius parameter 
%       s_idx:  s index.
%       d:      disturbances.
%               [d_ax; d_ay; d_w]^T 
% 
%   Returns:
%       dz_ds: 7x1 matrix. 
%              Error-state vector derivative in path-aligned frame:
%              [e_y'; e_phi'; v_x'; v_y'; w'; a_x'; a_y']^T
%   
%------------------------------------------------------------------    

    % ====== Unpack state ======
    e_y   = z(1);
    e_phi = z(2);
    v_x   = z(3);
    v_y   = z(4);
    w     = z(5);
    a_x   = z(6);
    a_y   = z(7);

    % ====== Unpack input ======
    alpha = u(1);
    j_x   = u(2);
    j_y   = u(3);

    % ====== Unpack disturbance (if not given, default 0) ======
    d_ax = d(1);
    d_ay = d(2);
    d_w  = d(3);

    % ====== Reference ======
    rho     = ref.rho;
    phi_dot = ref.phi_dot;

    % ====== Time-domain dz/dt ======
    % error dynamics
    e_y_p   = v_x * sin(e_phi) + v_y * cos(e_phi);
    e_phi_p = w - phi_dot;

    % velocities
    v_x_p   = a_x;
    v_y_p   = a_y;

    % accelerations + disturbance （這裡就是論文的 τ_x/m, τ_y/m, τ_z/J）
    a_x_p   = j_x   + d_ax;
    a_y_p   = j_y   + d_ay;
    w_p     = alpha + d_w;

    % ====== Compute s_dot ======
    v_s   = v_x * cos(e_phi) - v_y * sin(e_phi);
    s_dot = rho * v_s / (rho - e_y);   % 用你現在的定義

    eps_s = 1e-6;
    if abs(s_dot) < eps_s
        s_dot = sign(s_dot + eps_s) * eps_s;
    end

    dt_ds = 1 / s_dot;

    % ====== Convert to space-based dz/ds ======
    dz_ds = dt_ds * [ e_y_p;
                      e_phi_p;
                      v_x_p;
                      v_y_p;
                      w_p;
                      a_x_p;
                      a_y_p ];
end

