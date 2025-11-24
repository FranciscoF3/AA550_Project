function u = mpc_solver(z0, ref_vr, Ad, Bd, Cd, mpc)
% MPC_SOLVER  Quadprog-based MPC for space-based linearized model
%
%   z_{k+1} = Ad z_k + Bd u_k + Cd
%
%   States: z = [e_y; e_phi; v_x; v_y; w; a_x; a_y]
%   Inputs: u = [alpha; j_x; j_y]
%
%   Cost (same structure as CVX version):
%     J = sum_{k=1..Kh} [
%           W1*e_y(k)^2 + W2*e_phi(k)^2
%         + (v(k)-v_r)' W3 (v(k)-v_r)
%         + sigma(k)' W4 sigma(k)
%         + W5 * ||u_k - u_{k-1}||^2
%         ]
%         + W6 * (e_y(Kh+1)^2 + e_phi(Kh+1)^2)

    % ---------- Unpack MPC parameters ----------
    Kh = mpc.Kh;        % prediction horizon
    Ch = mpc.Ch;        % control horizon

    n = mpc.z_dim;      % 7
    m = mpc.u_dim;      % 3
    K = Kh;             % shorthand

    u_min = mpc.u_min;
    u_max = mpc.u_max;
    z_min = mpc.z_min;
    z_max = mpc.z_max;

    W1 = mpc.W1;        % ey
    W2 = mpc.W2;        % ephi
    W3 = mpc.W3;        % 3x3 on [vx; vy; w]
    W4 = mpc.W4;        % 2x2 on [sig_ey; sig_ephi]
    W5 = mpc.W5;        % Δu
    W6 = mpc.W6;        % terminal ey,ephi

    if isfield(mpc, 'u_prev')
        u_prev = mpc.u_prev;
    else
        u_prev = zeros(m,1);
    end

    v_r = ref_vr(:);    % [vx_ref; vy_ref; w_ref]

    % Decision variable x = [U_stack; s_stack]
    %   U_stack = [u_1; ...; u_K]        (m*K x 1)
    %   s_stack = [sigma_1; ...; sigma_K] (2*K x 1; sigma = [sig_ey; sig_ephi])
    Nu = m*K;
    Ns = 2*K;
    Nx = Nu + Ns;

    % ---------- Precompute A^k and sums for lifted model ----------
    % We want A_pows{k} = A^(k-1), k = 1..K+1
    % and S_pows{k} = sum_{i=0}^{k-1} A^i
    A_pows = cell(K+1,1);
    S_pows = cell(K+1,1);

    A_pows{1} = eye(n);      % A^0
    S_pows{1} = eye(n);      % I = sum_{i=0}^0 A^i

    for k = 2:K+1
        A_pows{k} = Ad * A_pows{k-1};          % A^(k-1)
        S_pows{k} = S_pows{k-1} + A_pows{k-1}; % sum_{i=0}^{k-1} A^i
    end

    % Helper: for step k (1..K+1), get F_k, G_k such that
    %   Z_k = G_k * U_stack + F_k
    function [F_k, G_k] = getFG(k)
        Ak = A_pows{k+1};    % A^k
        Sk = S_pows{k};      % sum_{i=0}^{k-1} A^i

        F_k = Ak * z0 + Sk * Cd;   % constant part

        G_k = zeros(n, Nu);
        % contribution from inputs u_1..u_min(k,K)
        for j = 1:min(k,K)
            Aj = A_pows{k-j+1};    % A^(k-j)
            colsU = (j-1)*m + (1:m);
            G_k(:, colsU) = Aj * Bd;
        end
    end

    % ---------- Build QP cost: 0.5 x'Hx + f'x ----------
    H = zeros(Nx, Nx);
    f = zeros(Nx, 1);

    % State weighting matrix for stage cost
    Q = zeros(n);
    Q(1,1)     = W1;           % ey
    Q(2,2)     = W2;           % ephi
    Q(3:5,3:5) = W3;           % vx, vy, w

    % 1) Stage costs for k = 1..K
    for k = 1:K
        [F_k, G_k] = getFG(k);

        % reference state [0;0;v_r]
        z_ref = zeros(n,1);
        z_ref(3:5) = v_r;

        E_k = F_k - z_ref;    % constant part in (Z_k - z_ref)

        % (G_k U + E_k)' Q (G_k U + E_k)
        QG  = Q * G_k;
        GQG = G_k' * QG;          % Nu x Nu
        GQE = G_k' * (Q * E_k);   % Nu x 1

        H(1:Nu, 1:Nu) = H(1:Nu, 1:Nu) + 2*GQG;
        f(1:Nu)       = f(1:Nu)       + 2*GQE;
    end

    % 2) Slack cost: sum_k sigma_k' W4 sigma_k
    for k = 1:K
        idx_s = Nu + (k-1)*2 + (1:2);  % indices of [sig_ey; sig_ephi] for stage k
        H(idx_s, idx_s) = H(idx_s, idx_s) + 2*W4;
    end

    % 3) Δu cost: W5 * sum_k ||u_k - u_{k-1}||^2
    for k = 1:K
        Ek = zeros(m, Nu);
        if k == 1
            Ek(:, 1:m) = eye(m);
            e_const = -u_prev;
        else
            cols_k  = (k-1)*m + (1:m);
            cols_km = (k-2)*m + (1:m);
            Ek(:, cols_k)  = eye(m);
            Ek(:, cols_km) = -eye(m);
            e_const = zeros(m,1);
        end

        EkT_Ek = Ek' * Ek;         % Nu x Nu
        EkT_ec = Ek' * e_const;    % Nu x 1

        H(1:Nu, 1:Nu) = H(1:Nu, 1:Nu) + 2*W5 * EkT_Ek;
        f(1:Nu)       = f(1:Nu)       + 2*W5 * EkT_ec;
    end

    % 4) Terminal cost: W6 * (ey_T^2 + ephi_T^2), T = K+1
    [F_T, G_T] = getFG(K+1);       % Z_{K+1}
    Q_T = zeros(n);
    Q_T(1,1) = W6;
    Q_T(2,2) = W6;

    E_T  = F_T;                    % error target = 0

    QG_T  = Q_T * G_T;
    GQG_T = G_T' * QG_T;
    GQE_T = G_T' * (Q_T * E_T);

    H(1:Nu, 1:Nu) = H(1:Nu, 1:Nu) + 2*GQG_T;
    f(1:Nu)       = f(1:Nu)       + 2*GQE_T;

    % ---------- Build inequality constraints Aineq x <= bineq ----------
    % 1) Input bounds u_min <= u_k <= u_max
    % 2) Soft corridor for ey,ephi (with slack)
    % 3) Hard bounds on z(3:7)
    % 4) Slack >= 0

    max_rows = (2*m + 4 + 2*(n-2))*K + Ns;
    Aineq = zeros(max_rows, Nx);
    bineq = zeros(max_rows, 1);
    row = 0;

    % 1) Input bounds
    for k = 1:K
        cols_u = (k-1)*m + (1:m);
        for i = 1:m
            % u_i(k) <= u_max(i)
            row = row+1;
            Aineq(row, cols_u(i)) =  1;
            bineq(row)            =  u_max(i);
            % u_i(k) >= u_min(i)  <=>  -u_i(k) <= -u_min(i)
            row = row+1;
            Aineq(row, cols_u(i)) = -1;
            bineq(row)            = -u_min(i);
        end
    end

    % 2) State bounds via Z_k(U) = G_k U + F_k
    for k = 1:K
        [F_k, G_k] = getFG(k);

        % slack indices at stage k
        idx_s1 = Nu + (k-1)*2 + 1;  % sig_ey
        idx_s2 = Nu + (k-1)*2 + 2;  % sig_ephi

        % ey soft: z_min(1) - sig1 <= e_y <= z_max(1) + sig1
        ey_rowU  = G_k(1,:);
        ey_const = F_k(1);

        % lower: e_y >= z_min(1) - sig1  => -ey_rowU*U - sig1 <= -z_min(1) + ey_const
        row = row+1;
        Aineq(row, 1:Nu)   = -ey_rowU;
        Aineq(row, idx_s1) = -1;
        bineq(row)         = -z_min(1) + ey_const;

        % upper: e_y <= z_max(1) + sig1  => ey_rowU*U - sig1 <= z_max(1) - ey_const
        row = row+1;
        Aineq(row, 1:Nu)   =  ey_rowU;
        Aineq(row, idx_s1) = -1;
        bineq(row)         =  z_max(1) - ey_const;

        % ephi soft: z_min(2) - sig2 <= e_phi <= z_max(2) + sig2
        ep_rowU  = G_k(2,:);
        ep_const = F_k(2);

        % lower
        row = row+1;
        Aineq(row, 1:Nu)   = -ep_rowU;
        Aineq(row, idx_s2) = -1;
        bineq(row)         = -z_min(2) + ep_const;

        % upper
        row = row+1;
        Aineq(row, 1:Nu)   =  ep_rowU;
        Aineq(row, idx_s2) = -1;
        bineq(row)         =  z_max(2) - ep_const;

        % hard bounds for z(3:7)
        for i = 3:n
            zi_rowU  = G_k(i,:);
            zi_const = F_k(i);

            % z_i >= z_min(i)  => -zi_rowU*U <= -z_min(i) + zi_const
            row = row+1;
            Aineq(row, 1:Nu) = -zi_rowU;
            bineq(row)       = -z_min(i) + zi_const;

            % z_i <= z_max(i)  => zi_rowU*U <= z_max(i) - zi_const
            row = row+1;
            Aineq(row, 1:Nu) =  zi_rowU;
            bineq(row)       =  z_max(i) - zi_const;
        end
    end

    % 4) Slack >= 0  =>  -sigma <= 0
    for i = 1:Ns
        row = row+1;
        Aineq(row, Nu+i) = -1;
        bineq(row)       = 0;
    end

    Aineq = Aineq(1:row,:);
    bineq = bineq(1:row);

    % ---------- Equality constraints: Aeq x = beq ----------
    % Control horizon: for k > Ch, U_k = U_Ch
    Aeq = [];
    beq = [];
    if Ch < K
        n_eq = (K-Ch)*m;
        Aeq = zeros(n_eq, Nx);
        beq = zeros(n_eq, 1);
        r = 0;
        for k = Ch+1:K
            cols_k  = (k-1)*m + (1:m);
            cols_ch = (Ch-1)*m + (1:m);
            for i = 1:m
                r = r+1;
                Aeq(r, cols_k(i))  =  1;
                Aeq(r, cols_ch(i)) = -1;
            end
        end
    end

    % ---------- Solve QP ----------
    H = (H + H')/2;  % enforce symmetry

    opts = optimoptions('quadprog', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point-convex');

    [x_opt, ~, exitflag] = quadprog(H, f, Aineq, bineq, Aeq, beq, [], [], [], opts);

    if exitflag <= 0 || any(isnan(x_opt))
        % fallback: use previous input if solver fails
        u = u_prev;
    else
        U_stack = x_opt(1:Nu);
        u = U_stack(1:m);
    end
end
