function [p_world, phi_world] = errorStateToWorld(z, ref_k)
%ERRORSTATETOWORLD  Convert error-state to world-frame pose (circle path).
%
%   Input:
%       z      : 7x1 error-state vector 
%                [e_y; e_phi; v_x; v_y; w; a_x; a_y]
%
%       ref_k  : struct containing reference info at step k:
%                ref_k.s       : arc-length at this step (unused here)
%                ref_k.rho     : circle radius (scalar)
%                ref_k.phi_r   : reference heading at this step (tangent)
%
%   Output:
%       p_world = [x; y]     : robot position in world frame
%       phi_world            : robot heading in world frame
%
%   Assumptions:
%       - Path is a circle of radius rho, centered at the origin (0,0).
%       - Heading phi_r is tangent direction (CCW motion).
%       - For a CCW circle: phi_r = theta + pi/2, so theta = phi_r - pi/2.
%

    % -------- unpack error state --------
    e_y   = z(1);
    e_phi = z(2);

    % -------- unpack reference --------
    R     = ref_k.rho;      % circle radius (scalar)
    phi_r = ref_k.phi_r;    % tangent heading at this s

    % -------- recover reference position on circle --------
    % theta: angle from x-axis to radius vector
    theta = phi_r - pi/2;   % for CCW motion

    x_r = R * cos(theta);
    y_r = R * sin(theta);

    % -------- path normal and world position --------
    % tangent: t_hat = [cos(phi_r); sin(phi_r)]
    % normal : n_hat = [-sin(phi_r); cos(phi_r)]
    n_hat = [-sin(phi_r);
              cos(phi_r)];

    p_world = [x_r; y_r] + e_y * n_hat;

    % -------- world heading --------
    phi_world = phi_r + e_phi;

end
