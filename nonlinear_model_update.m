function z_next = nonlin_model(u, z, delta_s, ref)
%/Update error state z(s) using the nonlinear space-based model.
% 
% Args:
%   u_cur:  3x1 matrix. 
%           Control input at current step:
%           [alpha; j_x; j_y]^T
%               alpha : yaw jerk (or angular acceleration command, depending on your model)
%               j_x   : longitudinal jerk
%               j_y   : lateral jerk
% 
%   z_cur:  7x1 matrix. 
%           Error-state vector in path-aligned frame:
%           [e_y; e_phi; v_x; a_x; v_y; a_y; w]^T
%               e_y   : lateral error
%               e_phi : heading error
%               v_x   : longitudinal velocity in robot frame
%               a_x   : longitudinal acceleration
%               v_y   : lateral velocity in robot frame
%               a_y   : lateral acceleration
%               w     : yaw rate
%          
%   delta_s: scalar.
%            Path step length Δs used for space-based integration.
% 
%   ref: struct.
%        Reference path information at current step, e.g.:
%           ref.phi_dot : reference heading rate (φ̇_r)
%           ref.tho     : path parameter (e.g. curvature-related term in ṡ definition)
% 
% Returns:
%   true_z_next: 7x1 matrix. 
%                Updated error-state at the next path step s + Δs:
%                [e_y; e_phi; v_x; a_x; v_y; a_y; w]^T
% /% 

    % ====== Derive params ======
    % State
    e_y   = z(1);
    e_phi = z(2);
    v_x   = z(3);
    a_x   = z(4);
    v_y   = z(5);
    a_y   = z(6);
    w     = z(7);

    % Control
    alpha = u(1);
    j_x   = u(2);
    j_y   = u(3);
    % ----------------------------

    % ======= Calculate f(z, u) ======
    % ------- Time-domain derivatives dz/dt -------
    e_y_p   = v_x * sin(e_phi) + v_y * cos(e_phi);
    e_phi_p = w - ref.phi_dot;
    v_x_p   = a_x;
    a_x_p   = j_x;
    v_y_p   = a_y;
    a_y_p   = j_y;
    w_p     = alpha; 
    
    % ------- ṡ and dt/ds -------
    v_s   = v_x * cos(e_phi) - v_y * sin(e_phi);           % robot vel along path
    s_dot = ref.tho * v_s / (ref.tho - e_y);               % ds/dt

    % Avoid division by zero / very small s_dot
    eps_s = 1e-6;
    if abs(s_dot) < eps_s
        s_dot = sign(s_dot) * eps_s;
    end

    s_dot_inv = 1 / s_dot;                                 % dt/ds
    
    % ------- Space-based dynamics dz/ds -------
    f = s_dot_inv *  [   e_y_p;
                         e_phi_p;
                         v_x_p;
                         a_x_p;
                         v_y_p;
                         a_y_p;
                         w_p;       ];

    % ====== Update z(s) ======
    z_next = z + f * delta_s;

end