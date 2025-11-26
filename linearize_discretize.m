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


function [A_r, B_r] = compute_jacobians(z_r, u_r, ref)
%   Compute_Jacobians

    % --- Log params ---
    % state
    ey    = z_r(1);
    ephi  = z_r(2);
    vx    = z_r(3);
    vy    = z_r(4);
    omega = z_r(5);
    ax    = z_r(6);
    ay    = z_r(7);

    % input
    alpha = u_r(1);
    jx    = u_r(2);
    jy    = u_r(3);

    % path radius
    rho_s = ref.rho;
    phi_dot_s = ref.phi_dot;

    % helper K
    K = vx * cos(ephi) - vy * sin(ephi);
    
    if abs(K) < 0.3
        K = 0.3;
    end

    nz = length(z_r);
    nu = length(u_r);

    A_r = zeros(nz);
    B_r = zeros(nz, nu);

    % f_1 = d e_y / ds
    % A_r(1, 1) = ((rho_s - ey) * (vx * sin(ephi) + vy * cos(ephi))^2) / (rho_s^2 * K^2) ...
               % - (vx * vy * (cos(ephi)^2 - sin(ephi)^2) + sin(ephi) * cos(ephi) * (vx^2 - vy^2)) / (rho_s * K^2);
    A_r(1, 1) = -(vx * sin(ephi) + vy * cos(ephi)) / (rho_s * K);
    % A_r(1, 2) = ((rho_s - ey) / rho_s) ...
    %            * (vx^2 + (vy^2 * (sin(ephi)^2 - cos(ephi)^2)) - (2 * vx * vy * sin(ephi) * cos(ephi))) / K^2;
    A_r(1, 2) = (rho_s - ey) * (vx^2 + vy^2) / (rho_s * K^2);

    A_r(1, 3) = (rho_s - ey) * (-vy) / (rho_s * K^2);
    A_r(1, 4) = (rho_s - ey) * vx / (rho_s * K^2);
    
    % f_2 = e_phi 
    A_r(2,1) = -(omega - phi_dot_s) / (rho_s * K);  % modefied
    A_r(2,2) = (omega - phi_dot_s) *(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi))/(rho_s*K^2);
    A_r(2, 3) = -(omega- phi_dot_s) * (rho_s - ey) * cos(ephi) / (rho_s * K^2);  % modefied
    A_r(2,4) = (omega - phi_dot_s)*(rho_s - ey)*sin(ephi)/(rho_s*K^2);
    A_r(2,5) = (rho_s - ey)/(rho_s*K);
    
    %%% f_3 = v_x %%
    A_r(3,1) = -ax/(rho_s*K);
    A_r(3,2) = ax*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi))/(rho_s*K^2);
    A_r(3,3) = -ax*(rho_s - ey)*cos(ephi)/(rho_s*K^2);
    A_r(3,4) = ax*(rho_s - ey)*sin(ephi)/(rho_s*K^2);
    A_r(3,6) = (rho_s - ey) / (rho_s * K);
    
    %%% f_4 = v_y %%
    A_r(4,1) = -ay/(rho_s*K);
    A_r(4,2) = ay*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2);
    A_r(4,3) = -ay*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(4,4) =  ay*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    A_r(4,7) = (rho_s - ey) / (rho_s*K);
    
    %%% f_5 = omega %%
    A_r(5,1) = -alpha/(rho_s*K);
    A_r(5,2) = alpha*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2);
    A_r(5,3) = -alpha*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(5,4) =  alpha*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    B_r(5,1) = (rho_s - ey) / (rho_s*K);
    
    %%% f_6 = ax %%
    A_r(6,1) = -jx/(rho_s*K);
    A_r(6,2) = jx*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2);
    A_r(6,3) = -jx*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(6,4) =  jx*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    B_r(6,2) =  (rho_s-ey) / (rho_s*K);

    %%% f_7 = a_y %%
    A_r(7,1) = -jy/(rho_s*K);
    B_r(7,3) = (rho_s - ey) / (rho_s * K);
    A_r(7,2) = jy*(rho_s - ey)*(vx*sin(ephi) + vy*cos(ephi)) / (rho_s*K^2); 
    A_r(7,3) = -jy*(rho_s - ey)*cos(ephi) / (rho_s*K^2);
    A_r(7,4) =  jy*(rho_s - ey)*sin(ephi) / (rho_s*K^2);
    
end
