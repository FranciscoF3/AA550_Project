function robot = robot_state_update(robot, z, ref, s_idx)
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
    
    % Log input
    ey    = z(1);
    ephi  = z(2);
    % vx    = z(3);
    % vy    = z(4);
    % omega = z(5);
    % ax    = z(6);
    % ay    = z(7);

    xr    = ref.x_r(s_idx);
    yr    = ref.y_r(s_idx);
    phi_r = ref.phi_r(s_idx);

    % Transfer error to global frame
    robot.x_cur = xr - ey * sin(phi_r);
    robot.y_cur = yr + ey * cos(phi_r);
    robot.phi_cur = phi_r + ephi;
    
end
