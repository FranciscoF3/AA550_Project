function test_jacobian_single()
    % test data (random)
    z = [0.05; 0.02; 0.8; 0.1; 0.1; 0.1; 0.1];
    u = [0.1;  0.1;  0.1];

    ref.rho     = 5.0;
    ref.phi_dot = 0.1;

    
    % test function in last doc
    addpath(fullfile('..'));


    [A_ana, B_ana] = compute_jacobians(z, u, ref);

    eps = 1e-6;
    A_num = numerical_jacobian_z(@(zz) f_continuous(zz, u, ref), z, eps);
    B_num = numerical_jacobian_u(@(uu) f_continuous(z, uu, ref), u, eps);

    A_err = A_ana - A_num;
    B_err = B_ana - B_num;

    fprintf('Max abs error in A: %.3e\n', max(abs(A_err(:))));
    fprintf('Max abs error in B: %.3e\n', max(abs(B_err(:))));

    disp('Row-wise max abs error in A:');
    disp(max(abs(A_err), [], 2));

    disp('Row-wise max abs error in B:');
    disp(max(abs(B_err), [], 2));

    % === Find max bias element in A ===
    [ maxA, idxA ] = max(abs(A_err(:)));
    [rowA, colA]   = ind2sub(size(A_err), idxA);
    
    fprintf('\nLargest A error at A(%d,%d): %.3e\n', rowA, colA, maxA);
    
    % === Find max bias element in B ===
    [ maxB, idxB ] = max(abs(B_err(:)));
    [rowB, colB]   = ind2sub(size(B_err), idxB);
    
    fprintf('Largest B error at B(%d,%d): %.3e\n', rowB, colB, maxB);

end
