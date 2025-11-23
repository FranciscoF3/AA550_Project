% AA 550- Nonlinear Optimal Control
% Project Report
% Authors: Kevin Ma, Luc Xavier Utheza, Francisco Flores



%%%// << Main MPC Loop >> //%%%

[ref] = reference_trajectory();
params = struct();
params.rho_s = ref.rho;

ds = ref.s(2) - ref.s(1);


%INitial robot position
robot.x_cur = ref.xr(1);
robot.y_cur = ref.yr(1);

% initial error state and input
z = zeros(7,1);
u = zeros(3,1);

steps = 300;
% Storage
traj_x = zeros(1,steps);
traj_y = zeros(1,steps);

for k = 1:steps

    %Step 1- Find the closest reference point
    ref_k = reference_sampling(ref, robot);
    
    z_r = z;
    u_r = u;

    %Step 1- linearize and discretize
    [Ad,Bd, Cd] = linearize_discretize(z_r, u_r, ds, params);

    %Step 3- compute control
    %%%%// << temp weights >> /// %%%%%
    % 3. MPC weights (simple LQR for now)
    Q = diag([20 20 1 1 0.1 0.1 0.1]);
    R = diag([0.1 0.1 0.1]);

    [K_lqr, ~,~] = dlqr(Ad, Bd, Q, R);
    u = -K_lqr*z;
    

    % step 4: update nonlinear error state
    z = nonlinear_step(z, u, ds, params);

    %Step 5- get robot state

    robot = robot_state_update(robot, z, ref_k);

    traj_x(k) = robot.x_cur;
    traj_y(k) = robot.y_cur;
end

figure;
hold on;
plot(ref.xr, ref.yr, '--k');
plot(traj_x, traj_y, 'r');
title('robot following circular path');



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
    dx = gradient(xr);
    dy = gradient(yr);
    phi_r = atan2(dy, dx);


    % Reference forward speed profile in s-plane
    v_r = 0.5 * ones(1, N);   % constant velocity with 0.5 m/s
    
    % Reference angular velocity
    phi_dot = v_r / R;

    % Save into structure
    ref.s     = s;
    ref.rho   = R;
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

function [Ad, Bd, Cd] = linearize_discretize(z_r, u_r, ds, params)
    
    [A_r, B_r, f_r] = compute_jacobians(z_r, u_r, params);
    
    %compute C_r
    C_r = f_r - A_r*z_r -B_r*u_r;
    
    nz = length(z_r);

    %%% Auxillary Matrix
    M_big = zeros(nz*2);

    % top left block A
    M_big(1:nz, 1:nz) = A_r;

    % top right block I
    M_big(1:nz, nz+1:2*nz) = eye(nz);

    M = expm(M_big *ds);

    M11 = M(1:nz,1:nz);
    M12 = M(1:nz, nz+1:2*nz);

    % Discrete matrices
    Ad = M11;   % A'=exp(A dt)
    Bd = M12 * B_r; %B'=M12*B
    Cd = M12 * C_r; % C'=M12*C
   
     
end

%%%
function [A_r, B_r, f_r] = compute_jacobians(z_r, u_r, params)

%Compute_Jacobians

    %state
    ey   = z_r(1);
    ephi = z_r(2);
    vx   = z_r(3);
    vy   = z_r(4);
    omega = z_r(5);
    ax   = z_r(6);
    ay   = z_r(7);

    %input
    alpha  = u_r(1);
    jx = u_r(2);
    jy = u_r(3);

    % path radius
    rho_s = params.rho_s;

    % helper K
    K = vx*cos(ephi)-vy*sin(ephi);

    nz = length(z_r);
    nu = length(u_r);

    A_r = zeros(nz);
    B_r = zeros(nz, nu);

    %%% f_1 = e_y %%
    A_r(1,1) = ((rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi))^2)/(rho_s^2*K^2) ...
               - (vx*vy*(cos(ephi)^2 - sin(ephi)^2) + sin(ephi)*cos(ephi)*(vx^2 - vy^2))/(rho_s*K^2);
    A_r(1,2) = ((rho_s - ey)*(vx^2 + vy^2)*(sin(ephi)^2 - cos(ephi)^2) ...
               - 2*vx*vy*sin(ephi)*cos(ephi))/K^2;
    A_r(1,3) = (rho_s - ey)/(rho_s*K^2);
    A_r(1,4) = -(vy)/(K^2);
    
    
    
    %%% f_2 = e_phi %%
    A_r(2,1) = -omega/(rho_s*K);
    A_r(2,2) = omega*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi))/(rho_s*K^2);
    A_r(2,3) = -omega*(rho_s - ey)*cos(ephi)/(rho_s*K^2);
    A_r(2,4) = omega*(rho_s - ey)*sin(ephi)/(rho_s*K^2);
    A_r(2,5) = (rho_s - ey)/(rho_s*K);
    
    %%% f_3 = v_x %%
    A_r(3,1) = -ax/(rho_s*K);
    A_r(3,2) = ax*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi))/(rho_s*K^2);
    A_r(3,3) = -ax*(rho_s - ey)*cos(ephi)/(rho_s*K^2);
    A_r(3,4) = ax*(rho_s - ey)*sin(ephi)/(rho_s*K^2);
    A_r(3,6) = (rho_s- ey)/ (rho_s*K);
    
    %%% f_4 = v_y %%
    A_r(4,1) = -ay/(rho_s*K);
    A_r(4,2) = ay*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2);
    A_r(4,3) = -ay*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(4,4) =  ay*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    A_r(4,7) = (rho_s - ey) / (rho_s*K);
    
    %%% f_5 = omega %%
    A_r(5,1) = -alpha/(rho_s*K);
    A_r(5,2) = alpha*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2);
    A_r(5,3) = -alpha*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(5,4) =  alpha*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    B_r(5,1) = (rho_s - ey) / (rho_s*K);
    
    %%% f_6 = ax %%
    A_r(6,1) = -jx/(rho_s*K);
    A_r(6,2) = jx*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2);
    A_r(6,3) = -jx*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(6,4) =  jx*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    B_r(6,2) =  (rho_s-ey) / (rho_s*K);

    %%% f_7 = a_y %%
    A_r(7,1) = -jy/(rho_s*K);
    B_r(7,3) = (rho_s - ey) / (rho_s*K);
    A_r(7,2) = jy*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2); 
    A_r(7,3) = -jy*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(7,4) =  jy*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    
 
    f_r = zeros(nz,1);
end


function robot = robot_state_update(robot, z, ref_k)
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

    ey = z(1);
    ephi = z(2);
    vx = z(3);
    vy = z(4);
    omega = z(5);
    ax = z(6);
    ay = z(7);

    xr = ref_k.xr;
    yr = ref_k.yr;
    phi_r = ref_k.phi_r;

    robot.x_cur = xr - ey*sin(phi_r);
    robot.y_cur = yr + ey*cos(phi_r);
    robot.phi_cur = phi_r + ephi;

    
    % --- continuous-time dynamics ---
    x_dot       = vx*cos(phi_r) - vy*sin(phi_r);
    y_dot       = vx*sin(phi_r) + vy*cos(phi_r);
    %phi_dot_dot = phi_ddot;

    vx_dot      = ax*cos(phi_r) - ay*sin(phi_r);
    vy_dot      = ax*sin(phi_r) + ay*cos(phi_r);

   % ax_dot      = jx;
    %ay_dot      = jy;


    robot.v_x_cur = x_dot;
    robot.v_y_cur = y_dot;
    robot.a_x_cur = vx_dot;
    robot.a_y_cur = vy_dot;
    robot.omega = omega;
end


function z_next = nonlinear_step(z, u, ds, params)
    ey   = z(1);
    ephi = z(2);
    vx   = z(3);
    vy   = z(4);
    omega = z(5);
    ax   = z(6);
    ay   = z(7);

    alpha = u(1);
    jx    = u(2);
    jy    = u(3);

    rho_s = params.rho_s;
    K = vx*cos(ephi) - vy*sin(ephi);

    % space-domain dynamics (very simple version)
    f1 = (vx*sin(ephi) + vy*cos(ephi)) / K;              % e_y'
    f2 = omega/rho_s - (vx*cos(ephi)-vy*sin(ephi))/K;    % e_phi'
    f3 = ax;                                             % v_x'
    f4 = ay;                                             % v_y'
    f5 = alpha;                                          % omega'
    f6 = jx;                                             % a_x'
    f7 = jy;                                             % a_y'

    f = [f1; f2; f3; f4; f5; f6; f7];

    z_next = z + ds * f;
end
