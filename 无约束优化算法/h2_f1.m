% 利用二分法，求解函数f(x) = x*sinx - x/6在[0,1]上的最小值，要求精度小于0.1，给出每一步的详细计算结果
f = @(x) x*sin(x) - x/6;
f_grad = @(x) sin(x) + x*cos(x) - 1/6;
a = 0, b = 1;
x = (a + b) / 2;
iter = 0;
while abs(b-a) > 0.1
    disp(['第',num2str(iter),'次迭代,当前x值为',num2str(x),',当前函数值为',num2str(f(x))])
    disp(['         当前a值为',num2str(a),'当前b值为',num2str(b)])
    if f_grad(x) > 0
        b = x;
    else
        a = x;
    end
    x = (a + b) / 2;
    iter = iter + 1;
end