% 利用黄金分割法，求解函数f(x) = x*sinx - x/6在[0,1]上的最小值，要求精度小于0.01，给出每一步的详细计算结果
f = @(x) x*sin(x) - x/6;
a = 0.0, b = 1.0;
p = (3 - sqrt(5)) / 2;
r = a + p*(b-a);
u = b - p*(b-a);
iter = 0;
while abs(b-a) > 0.1
    x = (a+b)/2;
    disp(['第',num2str(iter),'次迭代,当前区间中点值',num2str(x),',当前函数值为',num2str(f(x))])
    disp(['         当前区间为','[',num2str(a),',',num2str(b),']','，区间长度为',num2str(b-a)])
    if f(r) > f(u)
        a = r;
        r = u;
        u = b - p*(b-a);
    else
        b = u;
        u = r;
        r = a + p*(b-a);
    end  
    iter = iter + 1;
end
x = (a+b)/2;
disp(['第',num2str(iter),'次迭代,当前区间中点值',num2str(x),',当前函数值为',num2str(f(x))])
disp(['         当前区间为','[',num2str(a),',',num2str(b),']','，区间长度为',num2str(b-a)])