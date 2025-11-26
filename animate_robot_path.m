function animate_robot_path(X, Y, Phi, X_ref, Y_ref, t_log, robot, save_video, video_name)
%ANIMATE_ROBOT_PATH  播放 2D 上俯視的軌跡動畫
%
%   X, Y      : 1×N 機器人實際位置 (global frame)
%   Phi       : 1×N 機器人朝向 (rad, global frame)
%   X_ref,Y_ref : 1×N 參考軌跡
%   robot     : 結構，需包含 robot.L, robot.l（半長、半寬）
%   save_video: (optional) true/false 是否輸出 mp4
%   video_name: (optional) 檔名，例如 'omni_mpc_demo.mp4'

    if nargin < 8
        save_video = false;
    end
    if nargin < 9
        video_name = 'omni_mpc_demo.mp4';
    end

    N = numel(X);

    % 以 robot 幾何參數畫出車體外框（body frame 下）
    L = robot.L;   % 半長
    l = robot.l;   % 半寬

    body_shape = [...
        -L, -L,  L,  L;  % x 座標（車體座標）
        -l,  l,  l, -l]; % y 座標

    % ----- 建立圖形 -----
    figure; clf;
    hold on; grid on; axis equal;

    % 畫參考路徑
    plot(X_ref, Y_ref, 'k--', 'LineWidth', 1.2);

    % 設定顯示邊界（稍微多留一點 margin）
    margin = 1.0;
    xmin = min([X_ref(:); X(:)]) - margin;
    xmax = max([X_ref(:); X(:)]) + margin;
    ymin = min([Y_ref(:); Y(:)]) - margin;
    ymax = max([Y_ref(:); Y(:)]) + margin;
    axis([xmin xmax ymin ymax]);

    xlabel('X [m]');
    ylabel('Y [m]');
    title('Omni-directional Robot Path Tracking (MPC)');

    % 實際軌跡線
    h_traj = plot(NaN, NaN, 'b-', 'LineWidth', 1.5);

    % 車體 patch
    h_body = patch(NaN, NaN, [0.3 0.7 1.0]);  % 顏色可改

    % 前方朝向的小線（車頭）
    h_heading = plot(NaN, NaN, 'r-', 'LineWidth', 2);

    % 如要錄影，先準備 frame 陣列
    if save_video
        F(N) = struct('cdata',[], 'colormap',[]);
    end

    for k = 1:N
        % 更新軌跡線（從起點到目前）
        set(h_traj, 'XData', X(1:k), 'YData', Y(1:k));

        % 旋轉 + 平移車體外框
        R = [cos(Phi(k)) -sin(Phi(k));
             sin(Phi(k))  cos(Phi(k))];

        body_world = R * body_shape;
        bx = body_world(1,:) + X(k);
        by = body_world(2,:) + Y(k);
        set(h_body, 'XData', bx, 'YData', by);

        % 畫出車頭方向（例如車體中心往前畫一小段）
        head_len = L * 1.5;
        head_pt = [head_len; 0];   % 在 body frame 的一個點（前方）
        head_world = R * head_pt;
        hx = [X(k), X(k) + head_world(1)];
        hy = [Y(k), Y(k) + head_world(2)];
        set(h_heading, 'XData', hx, 'YData', hy);

        drawnow;

        if save_video
            F(k) = getframe(gcf);
        end

        if k > 1
            dt = t_log(k) - t_log(k-1);
            if dt > 0
                pause(dt);   % 真正依照物理時間播放
            end
        end
    end

    % ----- 輸出影片 -----
    if save_video
        v = VideoWriter(video_name, 'MPEG-4');
        v.FrameRate = 30;   % 播放速度
        open(v);
        writeVideo(v, F);
        close(v);
        fprintf('Video saved to %s\n', video_name);
    end
end
