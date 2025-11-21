% AA 550- Nonlinear Optimal Control
% Project Report
% Authors: Kevin Ma, Luc Xavier Utheza, Francisco Flores

params = struct();


function ref = reference_trajectory()
% Generate a reference circle path for MPC.
% 
% Returns:
%   ref: (struct) containing a discretized reference trajectory for a 
%   mobile robot. The path is defined as a circle with radius R in the 
%   global Cartesian frame. The trajectory is parameterized by the 
%   arc-length coordinate s, uniformly sampled with N points.
% 
%   The output struct 'ref' contains:
%       ref.s        - 1×N vector, arc-length samples along the trajectory.
%       ref.xr       - 1×N vector, x-coordinates in global frame.
%       ref.yr       - 1×N vector, y-coordinates in global frame.
%       ref.phi_r   - 1×N vector, reference heading (tangent direction).
%       ref.v_r      - 1×N vector, reference forward speed profile.
%       ref.phi_dot - 1×N vector, reference heading rate d(phi_r)/dt.

    % Path definition 
    R = 5.0;            % radius (meter)

    % Discretize path and project to s-space
    N = 1000;                   % sampling resolution of the reference table
    s = linspace(0, 2*pi*R, N); % projected mobile robot's position along the path

    % Coordinate transformation: s-plane to global Cartesian coordinates
    xr = R * cos(s / R);
    yr = R * sin(s / R);

    % Reference heading
    phi_r = atan2(diff([yr yr(end)]), diff([xr xr(end)]));

    % Reference forward speed profile in s-plane
    v_r = 0.5 * ones(1, N);   % constant velocity with 0.5 m/s
    
    % Reference angular velocity
    phi_dot = v_r / R;

    % Save into structure
    ref.s     = s;
    ref.tho   = R;
    ref.xr    = xr;
    ref.yr    = yr;
    ref.phi_r = phi_r;
    ref.v_r   = v_r;
    ref.phi_dot = phi_dot;
end

function target_ref = reference_sampling(ref, robot) 
%/Find the closest reference point on the desired path.
%
% Args: 
%   ref: reference path
%   x_cur, y_cur: (float) current position.
% 
% Returns:
%   target_ref: 
% /%

    % Compute distance between robot (x,y) and reference path
    dist = (ref.xr - robot.x_cur).^2 + (ref.yr - robot.y_cur).^2;

    % Find closest reference point
    [~, idx] = min(dist);

    target_ref.xr          = ref.xr(idx);
    target_ref.yr      = ref.yr(idx);
    target_ref.phi_r   = ref.phi_r(idx);
    target_ref.v_r     = ref.v_r(idx);
    target_ref.phi_dot = ref.phi_dot(idx);
end

function [Ad, Bd, Cd] = linearize_discretize(z_r, u_r, dt, params)

[A_r, B_r, f_r] = compute_jacobians(z_r, u_r, params);

C_r = f_r - A_r*z_r -B_r*u_r;

Ad = m

    % Extract blocks
    Ad = Md(1:nz, 1:nz);
    Bd = Md(1:nz, nz+1:nz+nu);
    Cd = Md(1:nz, nz+nu+1);
end

function [A_r, B_r, f_r] = compute_jacobians(z_r, u_r, params)

%Compute_Jacobians

K = v_x*cos(e_phi)-v_y*sin*(e_phi);

    A_r = 
    
    B_r =

    %%% f_1 = e_y %%
    A_r(1,1) = ((rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi))^2)/(rho_s^2*K^2) ...
               - (vx*vy*(cos(ephi)^2 - sin(ephi)^2) + sin(ephi)*cos(ephi)*(vx^2 - vy^2))/(rho_s*K^2);
    A_r(1,2) = ((rho_s - ey)*(vx^2 + vy^2)*(sin(ephi)^2 - cos(ephi)^2) ...
               - 2*vx*vy*sin(ephi)*cos(ephi))/K^2;
    A_r(1,3) = (rho_s - ey)/(rho_s*K^2);
    A_r(1,4) = -(vy)/(K^2);
    A_r(1,5) = -omega/(rho_s*K);
    
    
    %%% f_2 = e_phi %%
    A_r(2,1) = omega*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi))/(rho_s*K^2);
    A_r(2,2) = ((rho_s - ey)*(vx^2 + vy^2)*(sin(ephi)^2 - cos(ephi)^2) ...
                - 2*vx*vy*sin(ephi)*cos(ephi))/(K^2);
    A_r(2,3) = -omega*(rho_s - ey)*cos(ephi)/(rho_s*K^2);
    A_r(2,4) = omega*(rho_s - ey)*sin(ephi)/(rho_s*K^2);
    A_r(2,5) = (rho_s - ey)/(rho_s*K);
    
    %%% f_3 = v_x %%
    A_r(3,1) = -a_x/(rho_s*K);
    A_r(3,2) = a_x*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi))/(rho_s*K^2);
    A_r(3,3) = -a_x*(rho_s - ey)*cos(ephi)/(rho_s*K^2);
    A_r(3,4) = a_x*(rho_s - ey)*sin(ephi)/(rho_s*K^2);
    A_r(3,5) = (rho_s- e_y)/ (rho_s*K);
    
    %%% f_4 = v_y %%
    A_r(4,1) = -ay/(rho_s*K);
    A_r(4,2) = ay*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2);
    A_r(4,3) = -ay*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(4,4) =  ay*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    A_r(4,5) = (rho_s - ey) / (rho_s*K);
    
    %%% f_5 = omega %%
    A_r(5,1) = -alpha/(rho_s*K);
    A_r(5,2) = alpha*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2);
    A_r(5,3) = -alpha*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(5,4) =  alpha*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    B_r(5,1) = (rho_s - ey) / (rho_s*K);
    
    %%% f_6 = a_x %%
    A_r(6,1) = -jx/(rho_s*K);
    A_r(6,2) = jx*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2);
    A_r(6,3) = -jx*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(6,4) =  jx*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    A_r(6,5) =  (rho_s-ey) / (rho_s*K);

    %%% f_7 = a_y %%
    A_r(7,1) = -jy/(rho_s*K);
    A_r(7,2) = (rho_s - ey) / rho_sK;
    A_r(7,3) = jy*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2); 
    A_r(7,4) = -jy*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(7,5) =  jy*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    



end

function robot = robot_state_update(robot, u)
%/Update robot state through nonlinear dynamics model.
% 
% Args:
%   u: (arr) control input 3 x 10 matrix
%       j_x:   1x10 vector
%       j_y:   1x10 vector
%       alpha: 1x10 vector
% 
% Returns:
%   robot: (struct)-
% /%

    % --- continuous-time dynamics ---
    x_dot       = vx*cos(psi) - vy*sin(psi);
    y_dot       = vx*sin(psi) + vy*cos(psi);
    psi_dot_dot = psi_ddot;

    vx_dot      = ax + psi_dot*vy;
    vy_dot      = ay - psi_dot*vx;

    ax_dot      = jx;
    ay_dot      = jy;


    robot.x_cur = 
    robot.y_cur =
    robot.phi_cur = 
    robot.v_x_cur =
    robot.v_y_cur =
    robot.a_x_cur = 
    robot.a_y_cur = 
end
