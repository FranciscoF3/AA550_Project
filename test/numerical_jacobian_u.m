function J = numerical_jacobian_u(f_handle, u, eps)
% NUMERICAL_JACOBIAN_U  finite-difference Jacobian wrt input u
%
% f_handle: @(u) -> f(u)  (7x1)
% u       : current input (3x1)
% eps     : finite difference step

    nu = length(u);
    f0 = f_handle(u);      %#ok<NASGU>
    nf = length(f0);

    J = zeros(nf, nu);

    for i = 1:nu
        du      = zeros(nu,1);
        du(i)   = eps;

        f_plus  = f_handle(u + du);
        f_minus = f_handle(u - du);

        J(:,i) = (f_plus - f_minus) / (2*eps);
    end
end
