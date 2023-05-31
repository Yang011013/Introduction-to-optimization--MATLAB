% 共轭梯度法
% 定义函数
A = [3 0 1; 0 4 2; 1 2 3];
b = [3; 0; 1];
f = @(x) 0.5 * x' * A * x - b' * x;
grad_f = @(x) A * x - b;

x0 = [0; 0; 0];
d = -grad_f(x0);
alpha0 = -(d'*grad_f(x0)) / (d'*A*d);
x = x0 + alpha0 * d;
disp(['当前x:(',num2str(x(1)),',',num2str(x(2)),',',num2str(x(3)),')']);
disp(['函数值:', num2str(f(x))]);
while norm(grad_f(x)) > 0.01
    beta = (d'*A*grad_f(x)) / (d'*A*d);
    d = -grad_f(x) + beta*d;
    alpha = -(d'*grad_f(x)) / (d'*A*d);
    x = x + alpha*d;
    disp(['当前x:(',num2str(x(1)),',',num2str(x(2)),',',num2str(x(3)),')']);
    disp(['函数值:', num2str(f(x))]);
end