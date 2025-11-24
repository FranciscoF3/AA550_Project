function z_next = nonlinear_step(z, u, delta_s, ref)
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
    % Nonlinear continuous model
    dz_ds = f_continuous(z, u, ref);  

    % Update z(s)
    z_next = z + dz_ds * delta_s;   % Euler integration
   
end

