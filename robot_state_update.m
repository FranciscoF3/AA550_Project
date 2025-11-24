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
