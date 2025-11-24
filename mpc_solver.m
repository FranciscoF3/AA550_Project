function u = mpc_solver(z, zr, Ad, Bd, Cd, mpc)
%   MPC solver (CVX) for linearized space-based model
%
%   z_{k+1} = Ad z_k + Bd u_k + Cd
%
%   Args:
%       z   : (7x1 matrix). current error state (7x1)
%       z_r : (7x1 matrix). reference state used for linearization (7x1)
%       Ad,Bd,Cd : discrete-time linear model
%       mpc : struct with fields
%             mpc.Kh      - prediction horizon
%             mpc.Ch      - control horizon
%             mpc.z_dim   - state dimension (=7)
%             mpc.u_dim   - input dimension (=3)
%             mpc.W1,W2   - scalar weights for e_y, e_phi
%             mpc.W3      - 3x3 weight for [v_x; v_y; w]
%             mpc.W4      - scalar weight for slack sigma
%             mpc.u_min   - 3x1
%             mpc.u_max   - 3x1
%             mpc.e_y_min, mpc.e_y_max
%             mpc.e_phi_min, .e_phi_max

    % log params
    Kh = mpc.Kh;
    Ch = mpc.Ch;

    z_dim = mpc.z_dim; u_dim = mpc.u_dim; sigma_dim = mpc.sigma_dim;

    u_min = mpc.u_min; u_max = mpc.u_max;
    z_min = mpc.z_min; z_max = mpc.z_max;

    W1 = mpc.W1;
    W2 = mpc.W2;
    W3 = mpc.W3;
    W4 = mpc.W4;
    W5 = mpc.W5;
    W6 = mpc.W6;


    % Solve finite QP
    cvx_clear; 
    cvx_begin quiet
        % State z & Control input u
        variables Z(z_dim, Kh+1) U(u_dim, Kh) Sig(sigma_dim, Kh)

        % Initial condition
        Z(:, 1) == z;

        % Linearized prediction model
        for k = 1 : Kh
            Z(:, k+1) == Ad * Z(:, k) + Bd * U(:, k) + Cd;
        end

        % ------ Cost function ------
        J = 0;
        % reference velocity
        vr = zr(3:5);

        for k = 1 : Kh
            % errors
            e_y   = Z(1, k);
            e_phi = Z(2, k);

            % velocity state
            v_state = Z(3:5, k);

            J = J + W1 * e_y ^ 2 ...
                  + W2 * e_phi ^ 2 ...
                  + quad_form((v_state - vr), W3) ...
                  + quad_form(Sig(:, k), W4);  
        end

        minimize(J)

        % ------ Constraints ------
        subject to
            % input bounds & state corridor & slack
            for k = 1:Kh
                % input constraints
                U(:,k) >= u_min;
                U(:,k) <= u_max;

                % moving corridor (soften constraints)
                Z(1:5, k) >= (z_min(1:5) - Sig(:, k));
                Z(1:5, k) <= (z_max(1:5) + Sig(:, k));

                % slack nonnegativity
                Sig(:, k) >= 0;
                    
                % strict constraints
                Z(6:z_dim, k) >= z_min(6:z_dim); 
                Z(6:z_dim, k) <= z_max(6:z_dim);
            end 

            % --- control horizon constraint ---
            for k = Ch+1 : Kh
                U(:,k) == U(:,Ch);   % for k>=Ch, input = u_Ch
            end

    cvx_end

    u = U(:, 1);

end