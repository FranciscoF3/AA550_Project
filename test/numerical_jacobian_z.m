function J = numerical_jacobian_z(f_handle, z, eps)
% NUMERICAL_JACOBIAN_Z  finite-difference Jacobian wrt state z
%
% f_handle: @(z) -> f(z)  (7x1)
% z       : current state (7x1)
% eps     : finite difference step

    nz = length(z);
    f0 = f_handle(z);     %#ok<NASGU>   % 只是為了知道維度
    nf = length(f0);

    J = zeros(nf, nz);

    for i = 1:nz
        dz      = zeros(nz,1);
        dz(i)   = eps;

        f_plus  = f_handle(z + dz);
        f_minus = f_handle(z - dz);

        J(:,i) = (f_plus - f_minus) / (2*eps);
    end
end
