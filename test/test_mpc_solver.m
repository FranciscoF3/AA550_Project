clear; clc;

% 1. 讀取參數
[robot, mpc, sim] = init_params();

% 簡單一點，縮短 horizon 測試
mpc.Kh = 5;
mpc.Ch = 3;


% 2. 建一個「假的線性模型」，但尺寸正確
% z_{k+1} = z_k + B * u_k
z_dim = mpc.z_dim;
u_dim = mpc.u_dim;

Ad = eye(z_dim);               % 7x7，單位矩陣
Bd = zeros(z_dim, u_dim);      % 7x3，全零

% 讓 u(2) (j_x) 影響 v_x (state 3)
Bd(3,2) = 1;

Cd = zeros(z_dim, 1);          % 7x1

% 3. 設定初始狀態 & 參考
z   = zeros(z_dim, 1);   % 初始 error = 0
z(3) = -2;
mpc.W3 = 100 * eye(3);

z_r = zeros(z_dim, 1);

% 想要 v_x -> 1 [m/s]
z_r(3) = 1;

% 4. 呼叫 MPC
u = mpc_solver(z, z_r, Ad, Bd, Cd, mpc);

disp('Optimal control input u = ');
disp(u);
