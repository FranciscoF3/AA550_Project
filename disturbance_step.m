function d_next = disturbance_step(d_prev, dt, cfg)
%DISTURBANCE_STEP  Time-correlated random disturbance update.
%
%   d_next = disturbance_step(d_prev, dt, cfg)
%
%   d_prev : 3x1 previous disturbance [d_ax; d_ay; d_w]
%   dt     : time step [s]
%   cfg    : struct with fields
%            .tau_ax, .tau_ay, .tau_w   % time constants [s]
%            .sigma_ax, .sigma_ay, .sigma_w  % steady-state std dev
%
%   d_next : 3x1 updated disturbance

    % ---- default parameters ----
    if nargin < 3 || isempty(cfg)
        cfg.tau_ax   = 0.2;   % time constant [s]
        cfg.tau_ay   = 0.2;
        cfg.tau_w    = 0.2;
        cfg.sigma_ax = 0.10;  % steady-state std (調這個控制幅度)
        cfg.sigma_ay = 0.20;
        cfg.sigma_w  = 0.07;
    end

    if isempty(d_prev)
        d_prev = zeros(3,1);
    end

    % 一維一維更新（允許不同軸有不同特性）
    d_next = zeros(3,1);

    % x-direction
    alpha_ax = exp(-dt / cfg.tau_ax);
    q_ax     = sqrt(1 - alpha_ax^2) * cfg.sigma_ax;
    d_next(1) = alpha_ax * d_prev(1) + q_ax * randn;

    % y-direction
    alpha_ay = exp(-dt / cfg.tau_ay);
    q_ay     = sqrt(1 - alpha_ay^2) * cfg.sigma_ay;
    d_next(2) = alpha_ay * d_prev(2) + q_ay * randn;

    % yaw
    alpha_w = exp(-dt / cfg.tau_w);
    q_w     = sqrt(1 - alpha_w^2) * cfg.sigma_w;
    d_next(3) = alpha_w * d_prev(3) + q_w * randn;
end
