function [Ad, Bd, Cd] = linearize_discretize(z_r, u_r, ds, ref)
    % ------ Linearization -------
    [A_r, B_r] = compute_jacobians(z_r, u_r, ref);

    f_r = f_continuous(z_r, u_r, ref);
    
    C_r = f_r - A_r * z_r - B_r * u_r;

    % ------ Discretization ------
    nz = length(z_r);

    % Build M = exp([A I; 0 0] ds
    M_big = zeros(2*nz);
    M_big(1:nz,       1:nz      ) = A_r;
    M_big(1:nz,       nz+1:2*nz ) = eye(nz);
    % lower block is already zero

    M = expm(M_big * ds);
    M11 = M(1:nz, 1:nz);
    M12 = M(1:nz, nz+1:2*nz);

    Ad = M11;
    Bd = M12 * B_r;
    Cd = M12 * C_r;
   
end