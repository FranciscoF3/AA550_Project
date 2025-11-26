function z_ref = build_ref_with_obstacle(ref, mpc, sim)
%BUILD_REF_WITH_OBSTACLE
%   Pre-compute a global reference error-state profile z_ref(:,i)
%   over the entire path, including the effect of 1 static obstacle.
%
%   State definition:
%       z = [e_y;
%            e_phi;
%            v_x;
%            v_y;
%            w;
%            a_x;
%            a_y]
%
%   INPUT:
%       ref : struct of the reference path, with at least:
%             ref.s       : 1 x N, arc-length samples along the path
%             ref.x_r     : 1 x N, reference x (global)
%             ref.y_r     : 1 x N, reference y (global)
%             ref.v_r     : 1 x N, reference speed
%             ref.phi_dot : 1 x N, reference yaw rate
%
%       mpc : struct with obstacle info:
%             mpc.obs(1).x, mpc.obs(1).y    (global obstacle position)
%             OR mpc.obs(1).s0              (projected s-position)
%             mpc.obs(1).A                  (peak lateral offset in e_y)
%             mpc.obs(1).Delta              (influence half-width in s)
%
%       sim : struct with sim.ds (not strictly needed here)
%
%   OUTPUT:
%       z_ref : z_dim x N, reference error state along the path
%

    % ds currently not used, but keep if you want later
    ds    = sim.ds; %#ok<NASGU>
    N     = numel(ref.s);
    z_dim = mpc.z_dim;

    % Base reference: zero errors, track speed and yaw rate
    z_ref = zeros(z_dim, N);
    z_ref(3, :) = ref.v_r;      % v_x reference
    z_ref(5, :) = ref.phi_dot;  % w reference

    % No obstacle → just return the base reference
    if ~isfield(mpc, 'nObs') || mpc.nObs == 0
        return;
    end

    % For now assume exactly 1 static obstacle
    obs = mpc.obs(1);

    % -------------------------------------------------------------
    % 1) Find obstacle position along the path: s0 (in meters)
    % -------------------------------------------------------------
    if isfield(obs, 's0')
        % If you already know s0, just use it
        s0 = obs.s0;
    else
        % Project obstacle (x,y) onto path: find closest ref point
        x_o = obs.x;
        y_o = obs.y;

        dx = ref.x_r - x_o;
        dy = ref.y_r - y_o;
        d2 = dx.^2 + dy.^2;
        [~, idx_min] = min(d2);
        s0 = ref.s(idx_min);   % arc-length at closest point
    end

    % Parameters for lateral bump
    A     = obs.A;      % peak lateral offset [m]
    Delta = obs.Delta;  % influence half-width [m] along s

    if Delta <= 0
        return;
    end

    % -------------------------------------------------------------
    % 2) Build global cosine bump in e_y_ref(s)
    % -------------------------------------------------------------
    for i = 1:N
        ds_i = ref.s(i) - s0;   % s-distance from obstacle
        if abs(ds_i) <= Delta
            % Smooth cosine bump: 0 → A → 0 over [s0-Delta, s0+Delta]
            bump = 0.5 * A * (1 + cos(pi * ds_i / Delta));
            z_ref(1, i) = z_ref(1, i) + bump;   % e_y reference
        end
    end
end
