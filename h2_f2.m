% 利用Newton法，求解函数f(x) = x*sinx - x/6在[0,1]上的最小值，要求精度小于0.1，给出每一步的详细计算结果
f = @(x) x*sin(x) - x/6;
f_grad = @(x) sin(x) + x*cos(x) - 1/6;
f_grad_ = @(x) 2*cos(x) - x*sin(x);
next_x = 0.5;
current_x = 0.0;
iter = 0
while abs(next_x-current_x) > 0.1
    current_x = next_x;
    disp(['第',num2str(iter),'次迭代,当前x值为',num2str(current_x),',当前函数值为',num2str(f(current_x))])
    next_x = current_x - f_grad(current_x)/f_grad_(current_x);
    iter = iter + 1;
end
