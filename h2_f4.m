% 利用斐波那契数列法，求解函数f(x) = x*sinx - x/6在[0,1]上的最小值，要求精度小于0.1，给出每一步的详细计算结果
f = @(x) x*sin(x) - x/6;
% 通过计算需要得到的斐波那契数列长度为 n = 6
n = 6;
% 初始化斐波那契数列
F = [1 1];
for i = 3:n+1
    F(i) = F(i-1) + F(i-2);
end
disp([num2str(length(F))])
p = 1 - F(n)/F(n+1)
a = 0.0;
b = 1.0;
r = a + p*(b-a)
u = b - p*(b-a)
iter = 0
while abs(a-b) > 0.1
    x = (a+b)/2;
    disp(['第',num2str(iter+1),'次迭代,当前区间中点值',num2str(x),',当前函数值为',num2str(f(x))])
    disp(['         当前区间为','[',num2str(a),',',num2str(b),']','，区间长度为',num2str(b-a),'，当前p值为',num2str(p)])
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
    p = 1 - F(n-iter)/F(n+1-iter);
end
x = (a+b)/2;
disp(['第',num2str(iter+1),'次迭代,当前区间中点值',num2str(x),',当前函数值为',num2str(f(x))])
disp(['         当前区间为','[',num2str(a),',',num2str(b),']','，区间长度为',num2str(b-a),'，当前p值为',num2str(p)])
