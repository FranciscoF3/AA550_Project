function u = New_MPC_solver_QP(z0, zr_seq, Ad, Bd, Cd, mpc, ref_seq, s_idx)
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
%     + sigma_k' W4 sigma_k
%     + W5 * ||u_k - u_{k-1}||^2
%
%   Terminal cost:
%       W6 * (e_y(Kh+1)^2 + e_\phi(Kh+1)^2)
%
%   Slacks sigma_k soften corridor constraints on z(1:5).
%
%   INPUTS:
%       z0  : current state (7x1)
%       zr  : reference state for linearization (7x1)
%       Ad,Bd,Cd : discrete linear model (Kh fixed)
%       mpc : struct with fields:
%             Kh, Ch, z_dim, u_dim, sigma_dim
%             W1..W6, u_min/u_max, z_min/z_max
%             (optional) u_prev (3x1)
%
%   OUTPUT:
%       u   : optimal control input at current step (3x1)

    % ---------- MPC parameters ----------
    Kh   = mpc.Kh;         % prediction horizon
    Ch   = mpc.Ch;         % control horizon

    n    = mpc.z_dim;      % state dimension (7)
    m    = mpc.u_dim;      % input dimension (3)
    sdim = mpc.sigma_dim;  % number of slacks per stage (usually 5)

    % Constraints
    u_min = mpc.u_min;
    u_max = mpc.u_max;
    z_min = mpc.z_min;
    z_max = mpc.z_max;

    % Weights
    W1 = mpc.W1;   % ey
    W2 = mpc.W2;   % ephi
    W3 = mpc.W3;   % 3x3 on [vx; vy; w]
    W4 = mpc.W4;   % slack weight (scalar or sdim×sdim)
    W5 = mpc.W5;   % Δu penalty
    W6 = mpc.W6;   % terminal ey,ephi

    % Previous input for Δu term
    if isfield(mpc, 'u_prev')
        u_prev = mpc.u_prev;
    else
        u_prev = zeros(m,1);
    end

    % Reference velocity from zr
    % v_ref = zr(3:5);      % [v_x_ref; v_y_ref; w_ref]

    % % ---- Soft obstacle avoidance profile (ey_safe, Wobs) ----
    % % if mpc.obs is set，generate bump；or -> 0
    % if isfield(mpc, 'obs') && isfield(mpc, 'nObs') && mpc.nObs > 0
    %     [ey_safe, Wobs] = build_obstacle_bumps(ref_seq, mpc, s_idx);
    % else
    %     ey_safe = zeros(1, Kh);
    %     Wobs    = zeros(1, Kh);
    % end

    % ---------- Decision variable x = [Z_stack; U_stack; Sig_stack] ----------
    %   Z_stack = [z_1; ...; z_{Kh+1}]        (n*(Kh+1) x 1)
    %   U_stack = [u_1; ...; u_{Kh}]          (m*Kh x 1)
    %   Sig_stack = [sigma_1; ...; sigma_Kh]  (sdim*Kh x 1)
    NZ = n * (Kh+1);
    NU = m * Kh;
    NS = sdim * Kh;
    Nobs = 1 * Kh;
    Nx = NZ + NU + NS + Nobs;

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
    % z_ref       = zeros(n,1);
    % z_ref(3:5)  = v_ref;
    % Qzref       = Q * z_ref;

    N_ref = size(zr_seq, 2);   % total number of reference samples
    
    % Preallocate local reference for this horizon
    z_ref_seq = zeros(n, Kh+1);
    
    % Wrap-around slicing for circular path
    for k = 1:(Kh+1)
        idx = s_idx + (k-1);
        if idx > N_ref
            idx = idx - N_ref;   % wrap around for circle
        end
        z_ref_seq(:, k) = zr_seq(:, idx);
    end

    % z_ref_seq = zr_seq(:, s_idx: s_idx+Kh);

    % 1) Stage costs: state tracking + slack
    for k = 1:Kh

        iz = idxZ(k);
        is = idxS(k);

        % (z_k - z_ref)' Q (z_k - z_ref)
        % expands to z_k'Qz_k - 2 z_ref'Q z_k + const
        z_ref_k  = z_ref_seq(:, k);   
        Qzref  = Q * z_ref_k;
        H(iz,iz) = H(iz,iz) + 2*Q;
        f(iz)    = f(iz)    - 2*Qzref;

        % wobs_k = Wobs(k);
        % if wobs_k > 0
        %     idx_ey = iz(1);            
        %     ey_s   = ey_safe(k);       
        % 
        %     % Wobs(k) * (e_y - ey_s)^2
        %     H(idx_ey, idx_ey) = H(idx_ey, idx_ey) + 2*wobs_k;
        %     f(idx_ey)         = f(idx_ey)         - 2*wobs_k*ey_s;
        % end

        % sigma_k' W4 sigma_k
        H(is,is) = H(is,is) + 2*W4;

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
                % ||u_k - u_{k-1}||^2
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
    % 1) z_1 = z0
    % 2) z_{k+1} = Ad z_k + Bd u_k + Cd
    % 3) control horizon: for k > Ch, u_k = u_Ch

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
            % z_{k+1}(i) - Ad(i,:)*z_k - Bd(i,:)*u_k = Cd(i)
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
    % 2) Soft corridor for z(1:5) using slacks sigma
    % 3) Hard bounds on z(6:7)
    % 4) Slack >= 0

    max_rows = (2*m + 2*5 + 2*(n-5) + sdim + 1) * Kh;
    Aineq = zeros(max_rows, Nx);
    bineq = zeros(max_rows, 1);
    rI = 0;

    % 1) Input bounds
    for k = 1:Kh
        iu = idxU(k);
        for i = 1:m
            % u_i(k) <= u_max(i)
            rI = rI+1;
            Aineq(rI, iu(i)) =  1;
            bineq(rI)        =  u_max(i);
            % u_i(k) >= u_min(i) => -u_i(k) <= -u_min(i)
            rI = rI+1;
            Aineq(rI, iu(i)) = -1;
            bineq(rI)        = -u_min(i);
        end
    end

    % 2) Corridor with slacks on z(1:5), 3) hard bounds on z(6:7)
    for k = 1:Kh
        iz = idxZ(k);
        is = idxS(k);

        % softened constraints for z(1:5)
        for i = 1:5
            is_i = is(i);

            % lower: z_i >= z_min(i) - s_i  ->  -z_i - s_i <= -z_min(i)
            rI = rI+1;
            Aineq(rI, iz(i)) = -1;
            Aineq(rI, is_i)  = -1;
            bineq(rI)        = -z_min(i);

            % upper: z_i <= z_max(i) + s_i  ->   z_i - s_i <= z_max(i)
            rI = rI+1;
            Aineq(rI, iz(i)) =  1;
            Aineq(rI, is_i)  = -1;
            bineq(rI)        =  z_max(i);
        end

        % hard bounds for z(6:7)
        for i = 6:n
            % z_i >= z_min(i) -> -z_i <= -z_min(i)
            rI = rI+1;
            Aineq(rI, iz(i)) = -1;
            bineq(rI)        = -z_min(i);

            % z_i <= z_max(i)
            rI = rI+1;
            Aineq(rI, iz(i)) =  1;
            bineq(rI)        =  z_max(i);
        end


    end

    % 4) Slack >= 0 -> -sigma <= 0
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
    H = (H + H')/2;  % enforce symmetry

    opts = optimoptions('quadprog', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point-convex');

    % ====== DEBUG: check where complex comes from ======
    % if ~isreal(H)
    %     warning('H is complex. max imag = %g', max(abs(imag(H(:)))));
    % end
    % if ~isreal(f)
    %     warning('f is complex. max imag = %g', max(abs(imag(f(:)))));
    % end
    % if ~isreal(Aeq)
    %     warning('Aeq is complex. max imag = %g', max(abs(imag(Aeq(:)))));
    % end
    % if ~isreal(beq)
    %     warning('beq is complex. max imag = %g', max(abs(imag(beq(:)))));
    % end
    % if ~isreal(Aineq)
    %     warning('Aineq is complex. max imag = %g', max(abs(imag(Aineq(:)))));
    % end
    % if ~isreal(bineq)
    %     warning('bineq is complex. max imag = %g', max(abs(imag(bineq(:)))));
    % end


    [x_opt, ~, exitflag] = quadprog(H, f, Aineq, bineq, Aeq, beq, [], [], [], opts);

    % Fallback if solver fails
    if exitflag <= 0 || isempty(x_opt) || any(~isfinite(x_opt))
        warning('New_MPC_solver_QP: quadprog failed, using previous input.');
        u = u_prev;
    else
        % extract u_1 from U_stack
        U_stack = x_opt(NZ+1 : NZ+NU);
        u = U_stack(1:m);
    end
end
