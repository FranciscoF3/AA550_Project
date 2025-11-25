function u = mpc_solver(z, ref_vr, Ad, Bd, Cd, mpc)
%   MPC solver (CVX) for linearized space-based model
%
%   z      : current state (7x1)
%   ref_vr : desired velocity reference = [v_x_ref; v_y_ref; omega_ref]
%   Ad,Bd,Cd : discrete linearized error model
%   mpc    : parameter struct

    % Extract MPC parameters
    Kh = mpc.Kh;
    Ch = mpc.Ch;

    z_dim = mpc.z_dim; 
    u_dim = mpc.u_dim;

    % Slack: only 2 softened constraints: ey, ephi
    sigma_dim = 2;

    u_min = mpc.u_min; 
    u_max = mpc.u_max;
    z_min = mpc.z_min;
    z_max = mpc.z_max;

    W1 = mpc.W1;     % lateral error weight
    W2 = mpc.W2;     % heading error weight
    W3 = mpc.W3;     % velocity error penalty
    W4 = mpc.W4;     % slack penalty
    W5 = mpc.W5;     % Δu penalty
    W6 = mpc.W6;     % terminal penalty

    % Initial input (for Δu term)
    u_prev = mpc.u_prev;

    % ---- CVX MPC QP ----
    cvx_clear;
    cvx_begin quiet
        variables Z(z_dim, Kh+1) U(u_dim, Kh) Sig(sigma_dim, Kh)

        % Initial condition
        Z(:,1) == z;

        % System dynamics rollout
        for k = 1:Kh
            Z(:,k+1) == Ad*Z(:,k) + Bd*U(:,k) + Cd;
        end

        % Cost
        J = 0;
        v_r = ref_vr; % desired velocity [vx; vy; w]

        for k = 1:Kh
            e_y   = Z(1,k);
            e_phi = Z(2,k);

            v_state = Z(3:5,k);

            % base tracking costs
            J = J + ...
                W1 * e_y^2 + ...
                W2 * e_phi^2 + ...
                quad_form(v_state - v_r, W3) + ...
                quad_form(Sig(:,k), W4);

            % Δu rate penalty
            if k == 1
                du = U(:,k) - u_prev;
            else
                du = U(:,k) - U(:,k-1);
            end
            J = J + W5 * sum(du.^2);
        end

        % terminal cost (optional)
        J = J + W6 * (Z(1,Kh+1)^2 + Z(2,Kh+1)^2);

        minimize(J)

        % ---- Constraints ----
        subject to
            for k = 1:Kh

                % input limits
                U(:,k) >= u_min;
                U(:,k) <= u_max;

                % softened: ey, ephi
                Z(1,k) >= z_min(1) - Sig(1,k);
                Z(1,k) <= z_max(1) + Sig(1,k);

                Z(2,k) >= z_min(2) - Sig(2,k);
                Z(2,k) <= z_max(2) + Sig(2,k);

                % slack >= 0
                Sig(:,k) >= 0;

                % hard constraints on remaining states
                Z(3:7,k) >= z_min(3:7);
                Z(3:7,k) <= z_max(3:7);
            end

            % input freezing for k > Ch
            for k = Ch+1:Kh
                U(:,k) == U(:,Ch);
            end

    cvx_end

    % control action
    u = U(:,1);

end
