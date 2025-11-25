function u = New_MPC_solver_QP(z0, z_r, Ad, Bd, Cd, mpc)
% MPC_SOLVER  Quadprog-based MPC for space-based linearized model
%
%   z_{k+1} = Ad z_k + Bd u_k + Cd
%
%   States: z = [e_y; e_phi; v_x; v_y; w; a_x; a_y]
%   Inputs: u = [alpha; j_x; j_y]
%
%   Cost (per stage k = 1..Kh):
%       W1*e_y(k)^2 + W2*e_phi(k)^2
%     + (v(k)-v_ref)' W3 (v(k)-v_ref)
%     + Sig(:,k)' W4 Sig(:,k)
%     + W5 * ||u_k - u_{k-1}||^2
%   plus terminal cost W6 * (e_y(Kh+1)^2 + e_phi(Kh+1)^2).
%
%   Slacks Sig soften corridor on z(1:5).

    % ---------- Unpack MPC parameters ----------
    Kh   = mpc.Kh;
    Ch   = mpc.Ch;

    n    = mpc.z_dim;      % 7
    m    = mpc.u_dim;      % 3
    sdim = mpc.sigma_dim;  % usually 5

    u_min = mpc.u_min;
    u_max = mpc.u_max;
    z_min = mpc.z_min;
    z_max = mpc.z_max;

    W1 = mpc.W1;   % ey
    W2 = mpc.W2;   % ephi
    W3 = mpc.W3;   % 3x3 on [vx; vy; w]
    W4 = mpc.W4;   % scalar or sdim×sdim
    W5 = mpc.W5;   % Δu penalty
    W6 = mpc.W6;   % terminal ey,ephi

    % previous input for Δu (if provided)
    if isfield(mpc, 'u_prev')
        u_prev = mpc.u_prev;
    else
        u_prev = zeros(m,1);
    end

    % reference velocity from reference state z_r
    v_ref = z_r(3:5);   % [v_x_ref; v_y_ref; w_ref]

    % ---------- Decision variable x = [Z_stack; U_stack; Sig_stack] ----------
    %   Z_stack = [z_1; ...; z_{Kh+1}]      (n*(Kh+1) x 1)
    %   U_stack = [u_1; ...; u_{Kh}]        (m*Kh x 1)
    %   Sig_stack = [sigma_1; ...;sigma_Kh] (sdim*Kh x 1)
    NZ = n * (Kh+1);
    NU = m * Kh;
    NS = sdim * Kh;
    Nx = NZ + NU + NS;

    % index helpers into x
    idxZ = @(k) (k-1)*n + (1:n);                        % z_k
    idxU = @(k) NZ + (k-1)*m + (1:m);                   % u_k
    idxS = @(k) NZ + NU + (k-1)*sdim + (1:sdim);        % sigma_k

    % ---------- Build cost: 0.5 x'Hx + f'x ----------
    H = zeros(Nx, Nx);
    f = zeros(Nx, 1);

    % State weighting matrix for stage cost
    Q = zeros(n);
    Q(1,1)     = W1;          % ey
    Q(2,2)     = W2;          % ephi
    Q(3:5,3:5) = W3;          % vx, vy, w

    z_ref = zeros(n,1);
    z_ref(3:5) = v_ref;
    Qzref = Q * z_ref;

    % 1) Stage state + slack costs
    for k = 1:Kh
        iz = idxZ(k);
        is = idxS(k);

        % (z_k - z_ref)' Q (z_k - z_ref)
        H(iz,iz) = H(iz,iz) + 2*Q;
        f(iz)    = f(iz)    - 2*Qzref;

        % sigma_k' W4 sigma_k
        if isscalar(W4)
            H(is,is) = H(is,is) + 2*W4*eye(sdim);
        else
            H(is,is) = H(is,is) + 2*W4;
        end
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

    % 3) Terminal cost at z_{Kh+1}
    izT = idxZ(Kh+1);
    Q_T = zeros(n);
    Q_T(1,1) = W6;
    Q_T(2,2) = W6;
    H(izT,izT) = H(izT,izT) + 2*Q_T;

    % ---------- Equality constraints Aeq x = beq ----------
    % 1) z_1 = z0
    % 2) z_{k+1} = Ad z_k + Bd u_k + Cd
    % 3) control horizon: for k > Ch, U_k = U_Ch

    % base equalities: initial + dynamics
    nEq_base = n*(Kh+1);
    Aeq = zeros(nEq_base, Nx);
    beq = zeros(nEq_base, 1);
    row = 0;

    % 1) initial condition: Z(:,1) == z0
    iz1 = idxZ(1);
    for i = 1:n
        row = row+1;
        Aeq(row, iz1(i)) = 1;
        beq(row)         = z0(i);
    end

    % 2) dynamics: Z(:,k+1) == Ad*Z(:,k) + Bd*U(:,k) + Cd
    for k = 1:Kh
        izk  = idxZ(k);
        izk1 = idxZ(k+1);
        iuk  = idxU(k);

        for i = 1:n
            row = row+1;
            % z_{k+1}(i) - Ad(i,:)*z_k - Bd(i,:)*u_k = Cd(i)
            Aeq(row, izk1(i)) = 1;
            Aeq(row, izk)     = Aeq(row, izk) - Ad(i,:);
            Aeq(row, iuk)     = Aeq(row, iuk) - Bd(i,:);
            beq(row)          = Cd(i);
        end
    end

    % 3) Control horizon equalities: for k > Ch, U_k = U_Ch
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
    % 2) Soft corridor for z(1:5) with slacks Sig(:,k)
    % 3) Hard bounds on z(6:7)
    % 4) Slack >= 0

    max_rows = (2*m + 2*5 + 2*(n-5) + sdim) * Kh;
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

    % 2) Corridor with slacks for z(1:5), 3) hard bounds for z(6:7)
    for k = 1:Kh
        iz = idxZ(k);
        is = idxS(k);

        % softened z(1:5)
        for i = 1:5
            is_i = is(i);

            % lower: z_i >= z_min(i) - s_i -> -z_i - s_i <= -z_min(i)
            rI = rI+1;
            Aineq(rI, iz(i)) = -1;
            Aineq(rI, is_i)  = -1;
            bineq(rI)        = -z_min(i);

            % upper: z_i <= z_max(i) + s_i -> z_i - s_i <= z_max(i)
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
    H = (H + H')/2;  % ensure symmetry

    opts = optimoptions('quadprog', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point-convex');

    [x_opt, ~, exitflag] = quadprog(H, f, Aineq, bineq, Aeq, beq, [], [], [], opts);

    if exitflag <= 0 || isempty(x_opt) || any(~isfinite(x_opt))
        warning('mpc_solver (quadprog): solver failed, using previous input.');
        u = u_prev;
    else
        % extract u_1
        U_stack = x_opt(NZ+1 : NZ+NU);
        u = U_stack(1:m);
    end
end
