function animate_robot_path(X_log, Y_log, Phi_log, x_ref, y_ref, t_log, robot, mpc, save_video, video_name)
%ANIMATE_ROBOT_PATH  Play a 2D top-down animation of the robot trajectory.
%
%   X_log, Y_log: 1×N actual robot positions (global frame)
%   Phi_log     : 1×N robot heading (rad, global frame)
%   x_ref,y_ref : 1×N reference trajectory
%   robot       : struct containing robot.L and robot.l (half-length, half-width)
%   save_video  : (optional) true/false to export mp4
%   video_name  : (optional) filename such as 'omni_mpc_demo.mp4'

    if nargin < 9
        save_video = false;
    end
    if nargin < 10
        video_name = 'omni_mpc_demo.mp4';
    end

    N = numel(X_log);

    % Draw robot body outline (in body frame)
    L = robot.L;   % half-length
    l = robot.l;   % half-width

    body_shape = [...
        -L, -L,  L,  L;  % x-coordinates (body frame)
        -l,  l,  l, -l]; % y-coordinates

    % ----- Create figure -----
    figure; clf;
    hold on; grid on; axis equal;

    % Reference path
    p1 = plot(x_ref, y_ref, 'k--', 'LineWidth', 1.2);

    obs = mpc.obs(1);
    cx = obs.x;
    cy = obs.y;
    R  = obs.r;
    
    % plot obstacle center（with x）
    p2 = plot(cx, cy, 'rx', 'LineWidth', 2, 'MarkerSize', 10);
    
    % draw obstacle area
    theta = linspace(0, 2*pi, 200);
    x_circle = cx + R * cos(theta);
    y_circle = cy + R * sin(theta);
    p3 = plot(x_circle, y_circle, 'r-', 'LineWidth', 1.5);
    
    % fill with red
    patch(x_circle, y_circle, 'r', ...
          'FaceAlpha', 0.1, ...   
          'EdgeColor', 'none');
    % start / end
    plot(X_log(1),   Y_log(1),   'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot(X_log(end), Y_log(end), 'ms', 'MarkerSize', 10, 'MarkerFaceColor', 'm');

    % Display window (with extra margin)
    margin = 1.0;
    xmin = min([x_ref(:); X_log(:)]) - margin;
    xmax = max([x_ref(:); X_log(:)]) + margin;
    ymin = min([y_ref(:); Y_log(:)]) - margin;
    ymax = max([y_ref(:); Y_log(:)]) + margin;
    axis([xmin xmax ymin ymax]);
    
    legend([p1 p2, p3], {'Reference path', 'Obstacle center', 'Obstacle region'});

    xlabel('X [m]');
    ylabel('Y [m]');
    title('Omni-directional Robot Path Tracking (MPC)');

    % Actual trajectory line
    h_traj = plot(NaN, NaN, 'b-', 'LineWidth', 1.5);

    % Robot body patch
    h_body = patch(NaN, NaN, [0.3 0.7 1.0]); 

    % Heading direction indicator
    h_heading = plot(NaN, NaN, 'r-', 'LineWidth', 2);

    % Preallocate video frames
    if save_video
        F(N) = struct('cdata',[], 'colormap',[]);
    end

    for k = 1:N
        % Update trajectory (from start to current point)
        set(h_traj, 'XData', X_log(1:k), 'YData', Y_log(1:k));

        % Rotation + translation of robot body
        R = [cos(Phi_log(k)) -sin(Phi_log(k));
             sin(Phi_log(k))  cos(Phi_log(k))];

        body_world = R * body_shape;
        bx = body_world(1,:) + X_log(k);
        by = body_world(2,:) + Y_log(k);
        set(h_body, 'XData', bx, 'YData', by);

        % Heading line (extend from center forward)
        head_len = L * 2.2;
        head_pt = [head_len; 0];   % point in body frame
        head_world = R * head_pt;
        hx = [X_log(k), X_log(k) + head_world(1)];
        hy = [Y_log(k), Y_log(k) + head_world(2)];
        set(h_heading, 'XData', hx, 'YData', hy);

        drawnow;

        if save_video
            F(k) = getframe(gcf);
        end
        % Real-time playback based on log time
        if k > 1
            dt = t_log(k) - t_log(k-1);
            if dt > 0
                pause(dt);   
            end
        end
    end

    % ----- Export video -----
    if save_video
        v = VideoWriter(video_name, 'MPEG-4');
        v.FrameRate = 15;   
        open(v);
        writeVideo(v, F);
        close(v);
        fprintf('Video saved to %s\n', video_name);
    end
end
