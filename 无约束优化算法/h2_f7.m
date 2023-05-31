% 共轭梯度法
% 定义函数
A = [4 2; 2 2];
b = [-1; 1];
f = @(x) 0.5 * x' * A * x - b' * x;
grad_f = @(x) A * x - b;

% 定义共轭方向和初始点
d0 = [1; 0];
d1 = [-3/8; 3/4];
x0 = [0; 0];

% 计算步长 1
alpha = -(d0'*grad_f(x0)) / (d0'*A*d0);
x1 = x0 + alpha * d0;
disp(['当前x:(',num2str(x1(1)),',',num2str(x1(2)),')']);
disp(['函数值:', num2str(f(x1))]);

% 计算步长 2
alpha = -(d1'*grad_f(x1)) / (d1'*A*d1);
x2 = x1 + alpha * d1;
disp(['当前x:(',num2str(x2(1)),',',num2str(x2(2)),')']);
disp(['函数值:', num2str(f(x2))]);
   