function u = New_MPC_solver_QP(z0, zr, Ad, Bd, Cd, mpc, ref_seq)
% New_MPC_solver_QP
%   Quadprog-based MPC for space-based linearized model:
%
%       z_{k+1} = Ad z_k + Bd u_k + Cd
%
%   States: z = [e_y; e_\phi; v_x; v_y; w; a_x; a_y]
%   Inputs: u = [alpha; j_x; j_y]
%
%   Cost per stage k = 1..Kh:
%       W1*e_y(k)^2 + W2*e_\phi(k)^2
%     + (v(k)-v_ref)' W3 (v(k)-v_ref)
%     + sigma_corr(k)' W4_corridor sigma_corr(k)   (5 維 corridor slack)
%     + W4_obs * sigma_obs(k)^2                   (1 維 obstacle slack)
%     + W5 * ||u_k - u_{k-1||^2
%
%   Terminal cost:
%       W6 * (e_y(Kh+1)^2 + e_\phi(Kh+1)^2)
%
%   sigma_k : sdim×1, 其中
%       sigma_k(1:5) : corridor slack on z(1:5)
%       sigma_k(6)   : obstacle slack
%
%   mpc struct 除了原本欄位外，建議加：
%       mpc.W4_corridor (標量或 5x5)
%       mpc.W4_obs      (標量)
%   若未提供，則預設取 mpc.W4。
%

    % ---------- MPC parameters ----------
    Kh   = mpc.Kh;         % prediction horizon
    Ch   = mpc.Ch;         % control horizon

    n    = mpc.z_dim;      % state dimension (7)
    m    = mpc.u_dim;      % input dimension (3)
    sdim = mpc.sigma_dim;  % number of slacks per stage (>=6 for obstacle)

    if sdim < 6
        error('New_MPC_solver_QP: sigma_dim must be at least 6 (5 corridor + 1 obstacle).');
    end

    % Constraints
    u_min = mpc.u_min;
    u_max = mpc.u_max;
    z_min = mpc.z_min;
    z_max = mpc.z_max;

    % Weights
    W1 = mpc.W1;   % ey
    W2 = mpc.W2;   % ephi
    W3 = mpc.W3;   % 3x3 on [vx; vy; w]
    W5 = mpc.W5;   % Δu penalty
    W6 = mpc.W6;   % terminal ey,ephi

    % ---- Slack weights (不再混在同一個 W4 裡) ----
    if isfield(mpc, 'W4_corridor')
        W4_corr = mpc.W4_corridor;
    else
        % backward compatible: 用原本的 W4
        W4_corr = mpc.W4;
    end
    if isfield(mpc, 'W4_obs')
        W4_obs = mpc.W4_obs;
    else
        % 若沒給，預設比 corridor 再重一些
        if isscalar(mpc.W4)
            W4_obs = 10 * mpc.W4;
        else
            W4_obs = 10;   % 隨便給個合理值，建議你之後自己 tune
        end
    end

    % Previous input for Δu term
    if isfield(mpc, 'u_prev')
        u_prev = mpc.u_prev;
    else
        u_prev = zeros(m,1);
    end

    % ------------ obstacle params ----------------
    obs_x = mpc.obs.cx;
    obs_y = mpc.obs.cy;
    obs_r = mpc.obs.R_safe;
    nx    = mpc.obs.nx;
    ny    = mpc.obs.ny;

    % reference path (world frame)
    xr_vec   = ref_seq.x_r(:);    % Kh×1
    yr_vec   = ref_seq.y_r(:);
    phi_vec  = ref_seq.phi_r(:);

    % Reference velocity from zr
    v_ref = zr(3:5);      % [v_x_ref; v_y_ref; w_ref]

    % ---------- Decision variable x = [Z_stack; U_stack; Sig_stack] ----------
    NZ = n * (Kh+1);
    NU = m * Kh;
    NS = sdim * Kh;
    Nx = NZ + NU + NS;

    % index helpers into x
    idxZ = @(k) (k-1)*n + (1:n);                 % indices for z_k
    idxU = @(k) NZ + (k-1)*m + (1:m);            % indices for u_k
    idxS = @(k) NZ + NU + (k-1)*sdim + (1:sdim); % indices for sigma_k

    % ---------- Build cost: 0.5 x'Hx + f'x ----------
    H = zeros(Nx, Nx);
    f = zeros(Nx, 1);

    % State tracking weight Q
    Q = zeros(n);
    Q(1,1)     = W1;       % ey
    Q(2,2)     = W2;       % ephi
    Q(3:5,3:5) = W3;       % [vx, vy, w]

    % reference full-state vector used in (z_k - z_ref)
    z_ref      = zeros(n,1);
    z_ref(3:5) = v_ref;
    Qzref      = Q * z_ref;

    % 1) Stage costs: state tracking + slack (拆成 corridor / obstacle)
    for k = 1:Kh
        iz = idxZ(k);
        is = idxS(k);

        % (z_k - z_ref)' Q (z_k - z_ref)
        H(iz,iz) = H(iz,iz) + 2*Q;
        f(iz)    = f(iz)    - 2*Qzref;

        % ---- corridor slack sigma(1:5) ----
        is_corr = is(1:5);
        if isscalar(W4_corr)
            H(is_corr, is_corr) = H(is_corr, is_corr) + 2*W4_corr*eye(5);
        else
            % 如果你給的是 6x6 或更大，只取前 5x5 給 corridor 用
            Wc = W4_corr;
            if all(size(Wc) == [sdim sdim])
                Wc = Wc(1:5, 1:5);
            elseif ~all(size(Wc) == [5 5])
                error('W4_corridor must be scalar, 5x5, or sdim x sdim.');
            end
            H(is_corr, is_corr) = H(is_corr, is_corr) + 2*Wc;
        end

        % ---- obstacle slack sigma(6) ----
        s_obs = is(6);
        % W4_obs 一般是 scalar
        H(s_obs, s_obs) = H(s_obs, s_obs) + 2*W4_obs;
    end

    % 2) Δu penalty: W5 * sum_k ||u_k - u_{k-1}||^2
    if W5 > 0
        for k = 1:Kh
            iu = idxU(k);
            if k == 1
                % ||u_1 - u_prev||^2
                H(iu,iu) = H(iu,iu) + 2*W5*eye(m);
                f(iu)    = f(iu)    - 2*W5*u_prev;
            else
                iu_prev = idxU(k-1);
                H(iu,iu)           = H(iu,iu)           + 2*W5*eye(m);
                H(iu_prev,iu_prev) = H(iu_prev,iu_prev) + 2*W5*eye(m);
                H(iu,iu_prev)      = H(iu,iu_prev)      - 2*W5*eye(m);
                H(iu_prev,iu)      = H(iu_prev,iu)      - 2*W5*eye(m);
            end
        end
    end

    % 3) Terminal cost: W6 * (ey_T^2 + ephi_T^2) at z_{Kh+1}
    izT = idxZ(Kh+1);
    Q_T = zeros(n);
    Q_T(1,1) = W6;
    Q_T(2,2) = W6;
    H(izT,izT) = H(izT,izT) + 2*Q_T;

    % ---------- Equality constraints Aeq x = beq ----------
    nEq_base = n*(Kh+1);
    Aeq = zeros(nEq_base, Nx);
    beq = zeros(nEq_base, 1);
    row = 0;

    % 1) Initial condition: z_1 = z0
    iz1 = idxZ(1);
    for i = 1:n
        row = row+1;
        Aeq(row, iz1(i)) = 1;
        beq(row)         = z0(i);
    end

    % 2) Dynamics: z_{k+1} = Ad z_k + Bd u_k + Cd
    for k = 1:Kh
        izk  = idxZ(k);
        izk1 = idxZ(k+1);
        iu   = idxU(k);

        for i = 1:n
            row = row+1;
            Aeq(row, izk1(i)) = 1;
            Aeq(row, izk)     = Aeq(row, izk) - Ad(i,:);
            Aeq(row, iu)      = Aeq(row, iu)  - Bd(i,:);
            beq(row)          = Cd(i);
        end
    end

    % 3) Control horizon: for k > Ch, u_k = u_Ch
    Aeq_ch = [];
    beq_ch = [];
    if Ch < Kh
        nEq_ch = (Kh - Ch)*m;
        Aeq_ch = zeros(nEq_ch, Nx);
        beq_ch = zeros(nEq_ch, 1);
        r = 0;
        iu_ch = idxU(Ch);
        for k = Ch+1:Kh
            iu_k = idxU(k);
            for i = 1:m
                r = r+1;
                Aeq_ch(r, iu_k(i))  =  1;
                Aeq_ch(r, iu_ch(i)) = -1;
            end
        end
    end

    Aeq = [Aeq; Aeq_ch];
    beq = [beq; beq_ch];

    % ---------- Inequality constraints Aineq x <= bineq ----------
    % 1) Input bounds u_min <= u_k <= u_max
    % 2) Soft corridor for z(1:5) using slacks sigma(1:5)
    % 3) Hard bounds on z(6:7)
    % 4) Slack >= 0
    % 5) Obstacle avoidance with slack sigma(6)

    max_rows = (2*m + 2*5 + 2*(n-5) + 1 + sdim) * Kh;
    Aineq = zeros(max_rows, Nx);
    bineq = zeros(max_rows, 1);
    rI = 0;

    % 1) Input bounds
    for k = 1:Kh
        iu = idxU(k);
        for i = 1:m
            rI = rI+1;
            Aineq(rI, iu(i)) =  1;
            bineq(rI)        =  u_max(i);

            rI = rI+1;
            Aineq(rI, iu(i)) = -1;
            bineq(rI)        = -u_min(i);
        end
    end

    % 2) Corridor + 3) hard bounds + 5) obstacle
    for k = 1:Kh
        iz = idxZ(k);
        is = idxS(k);

        % corridor slacks on z(1:5)
        for i = 1:5
            s_idx = is(i);

            % lower: z_i >= z_min(i) - s_i  -> -z_i - s_i <= -z_min(i)
            rI = rI+1;
            Aineq(rI, iz(i)) = -1;
            Aineq(rI, s_idx) = -1;
            bineq(rI)        = -z_min(i);

            % upper: z_i <= z_max(i) + s_i  ->  z_i - s_i <= z_max(i)
            rI = rI+1;
            Aineq(rI, iz(i)) =  1;
            Aineq(rI, s_idx) = -1;
            bineq(rI)        =  z_max(i);
        end

        % hard bounds for z(6:7)
        for i = 6:n
            rI = rI+1;
            Aineq(rI, iz(i)) = -1;
            bineq(rI)        = -z_min(i);

            rI = rI+1;
            Aineq(rI, iz(i)) =  1;
            bineq(rI)        =  z_max(i);
        end

        % obstacle avoidance with slack sigma_obs = sigma(6)
        xr_k  = xr_vec(k);
        yr_k  = yr_vec(k);
        phi_k = phi_vec(k);

        c0 = nx*(xr_k - obs_x) + ny*(yr_k - obs_y);
        c1 = -nx*sin(phi_k) + ny*cos(phi_k);

        % 避免係數太小造成數值病態
        if abs(c1) < 1e-3
            c1 = sign(c1 + 1e-6) * 1e-3;
        end

        s_obs = is(6);  % 第六維 slack

        % c0 + c1*e_y >= obs_r - sigma_obs
        % -> -c1*e_y - sigma_obs <= c0 - obs_r
        rI = rI+1;
        Aineq(rI, iz(1)) = -c1;      % e_y
        Aineq(rI, s_obs) = -1;       % -sigma_obs
        bineq(rI)        = c0 - obs_r;
    end

    % 4) Slack >= 0 -> -sigma <= 0  (所有 slack 都非負)
    for k = 1:Kh
        is = idxS(k);
        for i = 1:sdim
            rI = rI+1;
            Aineq(rI, is(i)) = -1;
            bineq(rI)        = 0;
        end
    end

    Aineq = Aineq(1:rI,:);
    bineq = bineq(1:rI);

    % ---------- Solve QP ----------
    % 強制全部都是實數，避免 tiny imaginary 傳下去害 X_log/Y_log 變 complex
    H     = real((H + H')/2);
    f     = real(f);
    Aeq   = real(Aeq);
    beq   = real(beq);
    Aineq = real(Aineq);
    bineq = real(bineq);

    opts = optimoptions('quadprog', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point-convex');

    [x_opt, ~, exitflag] = quadprog(H, f, Aineq, bineq, Aeq, beq, [], [], [], opts);

    % Fallback if solver fails
    if exitflag <= 0 || isempty(x_opt) || any(~isfinite(x_opt))
        warning('New_MPC_solver_QP: quadprog failed, using previous input.');
        u = u_prev;
    else
        U_stack = x_opt(NZ+1 : NZ+NU);
        u = U_stack(1:m);
    end
end
