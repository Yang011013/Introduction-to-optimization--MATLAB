[TOC]

##  作业一 线性代数基础

### Q1

提交一个函数f1，传参为m、n (m、n为正整数)，生成一个m*n的矩阵，矩阵的所有值均为1

```matlab
% 传入参数为m,n，生成m*n的矩阵，矩阵所有元素为1
function f1(m,n)
	matrix = ones(m,n)
end
```

### Q2

提交一个函数f2，传参为m、n (m、n为正整数)，生成一个m*n的矩阵，矩阵的所有值均为3

```matlab
% 传参为m、n (m、n为正整数)，生成一个m*n的矩阵，矩阵的所有值均为3
function f2(m,n)
	matrix = ones(m,n) * 3 
end
```

### Q3

提交一个函数f3，传参为一个3*3的矩阵，返回该矩阵的逆矩阵，如果无逆矩阵，返回-1

```matlab
% 传参为一个3*3的矩阵，返回该矩阵的逆矩阵，如果无逆矩阵，返回-1
function flag = f3(A)
    if det(A) == 0 % 行列式等于0，矩阵不可逆
        flag = -1
    else
        flag = inv(A)
    end
end
```

### Q4

提交一个函数f4，传参为一个5*5的矩阵，判断该矩阵是否为对称矩阵，如果是则返回1，否则返回0

```matlab
% 传参为一个5*5的矩阵，判断该矩阵是否为对称矩阵，如果是则返回1，否则返回0
function flag = f4(A)
    if A' == A
        flag = 1
    else
        flag = 0      
end
```

### Q5

提交一个函数f5，传参为3个长度为3的列向量，判断三个列向量是否线性相关，如果是则返回1，否则返回0

```matlab
% 传参为3个长度为3的列向量，判断三个列向量是否线性相关，如果是则返回1，否则返回0
function flag = f5(a1,a2,a3)
    A = [a1; a2; a3]
    if det(A) == 0
        flag = 1
    else
        flag = 0
end
```

### Q6

提交一个函数f6，传参为一个10*10的矩阵、m、n (m、n为正整数且小于等于10)，返回该矩阵第m、n行交换后的矩阵的转置，要求用初等行变换的方法完成

```matlab
% 传参为一个10*10的矩阵、m、n (m、n为正整数且小于等于10)
% 返回该矩阵第m、n行交换后的矩阵的转置，要求用初等行变换的方法完成
function A = f6(A, m, n)
    A([m,n],:) = A([n,m],:)
    A = A'
end
```

### Q7

提交一个函数f7，传参为一个10*10的矩阵、m、k (m正整数且小于等于10)，返回该矩阵第m行乘以k后的矩阵的转置，要求用初等行变换的方法完成

```matlab
% 传参为一个10*10的矩阵、m、k (m正整数且小于等于10)
% 返回该矩阵第m行乘以k后的矩阵的转置，要求用初等行变换的方法完成
function A = f7(A,m,k)
    A(m,:) = A(m,:) * k
    A = A'
end
```

### Q8

提交一个函数f8，传参为一个10*10的矩阵、m、n、k (m、n为正整数且小于等于10)，返回该矩阵第m行乘以k加到第n行后的矩阵的转置，要求用初等行变换的方法完成

```matlab
% 传参为一个10*10的矩阵、m、n、k (m、n为正整数且小于等于10)
% 返回该矩阵第m行乘以k加到第n行后的矩阵的转置，要求用初等行变换的方法完成
function A = f8(A,m,n,k)
    A(n,:) = A(n,:) +  A(m,:) * k
    A = A'
end
```

### Q9

提交一个函数f9，传参为三个实数a、b、c (a>0)，绘制一元二次函数y=a*x*x+b*x+c的函数曲线示意图，并在图中显示出最小值的横纵坐标

```matlab
% 传参为三个实数a、b、c (a>0)
% 绘制一元二次函数y=a*x*x+b*x+c的函数曲线示意图，并在图中显示出最小值的横纵坐标
function f9(a,b,c)
    if a <= 0
        error('a必须大于0')
    end
    f = @(x) a*x.^2 + b.*x + c;
    x = -5:0.1:5;
    plot(x,f(x))
    hold on
    min_x = -b / (2 * a)
    min_f = f(min_x)
    plot(min_x,min_f,'r*')
    hold on
    text(min_x,min_f,'最小值坐标')
end
```

## 作业二 无约束优化问题

### Q1

利用二分法，求解函数$f(x) = x*sinx - x/6$在$[0,1]$上的最小值，要求精度小于0. 1，给出每一步的详细计算结果；

```matlab
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
```

### Q2

利用Newton法，求解函数$f(x) = x*sinx - x/6$在$[0,1]$上的最小值，要求精度小于0. 1，给出每一步的详细计算结果；

```matlab
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
```

### Q3

利用黄金分割法，求解函数$f(x) = x*sinx - x/6$在$[0,1]$上的最小值，要求精度小于0. 1，给出每一步的详细计算结果；

```matlab
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
```

### Q4

利用斐波那契数列法，求解函数$f(x) = x*sinx - x/6$在$[0,1]$上的最小值，要求精度小于0. 1，给出每一步的详细计算结果；

计算n，即迭代次数$(1+2*\epsilon)/F_{n+1}\leq n; F_{-1}=0,F_0=1,F_1=1,F_2=2,F_3=3,F_4=5 ...$,

```matlab
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
```

以上都是一元的函数。下面是多元的无约束优化问题。

### Q5

利用最速下降法，求解二元函数$f(x,y) = x^4+2*y^4在x\in [-4,4],y\in [-4,4]$上的最小值，起始点为$x=4,y=4$，要求精度小于0.1，给出每一步的计算结果。

函数工具：`fminbnd`

```matlab
% 利用最速下降法，求解二元函数f(x,y) = x^4+2*y^4在x\in [-4,4],y\in [-4,4]上的最小值
% 起始点为x=4,y=4，要求精度小于0.1，给出每一步的计算结果。
f = @(x,y) x^4 + 2*y^4;
grad_f = @(x,y) [4*x^3, 8*y^3];

start_point = [4,4];
precision = 0.1;
current_point = start_point;
iter = 0;
disp(['迭代次数 0：']);
disp(['当前点：（',num2str(current_point(1)),',',num2str(current_point(2)),')']);
disp(['当前函数值：',num2str(f(current_point(1),current_point(2)))])
% 迭代更新
while true
    gradient = grad_f(current_point(1),current_point(2)); 
    % 精确搜索步长，这里是对于alpha的医院函数，可以使用一元函数优化的算法，例如牛顿法，黄金分割法，二分法等。此处用的MATLAB工具箱的函数fminbnd
    % 搜索区间[-100,100]最小值
    alpha = fminbnd(@(alpha)f(current_point(1)-alpha*gradient(1),current_point(2)-alpha*gradient(2)),-100,100); 
    next_point = current_point - alpha * gradient;
    iter = iter + 1;
    disp(' ')
    disp(['迭代次数',num2str(iter),':']);
    disp(['迭代步长alpha:',num2str(alpha)]);
    disp(['当前点：（',num2str(next_point(1)),',',num2str(next_point(2)),')']);
    disp(['当前函数值：',num2str(f(next_point(1),next_point(2)))]);
    if norm(next_point - current_point) < precision
        break
    end
    current_point = next_point;
end
disp(' ')
disp(['迭代完成']);
disp(['最小值对应的点：（',num2str(next_point(1)),',',num2str(next_point(2)),')']);
disp(['函数最小值：',num2str(f(next_point(1),next_point(2)))]);
```

### Q6

利用Newton法，求解二元函数$f(x,y) = x^4+2*y^4在x\in [-4,4],y\in [-4,4]$上的最小值，起始点为$x=4,y=4$，要求精度小于0.1，给出每一步的计算结果。

```matlab
% 利用Newton法，求解二元函数f(x,y) = x^4+2*y^4在x\in [-4,4],y\in [-4,4]上的最小值
% 起始点为x=4,y=4，要求精度小于0.1，给出每一步的计算结果。
% 定义函数
f = @(x,y) x^4 + 2*y^4;
% 定义梯度向量
grad_f = @(x, y) [4*x^3; 8*y^3];
% 定义海森矩阵
hessian = @(x, y) [12*x^2, 0; 0, 24*y^2];

start_point = [4,4];
precision = 0.1;
current_point = start_point;
iter = 0;
disp(['迭代次数 0：']);
disp(['当前点：（',num2str(current_point(1)),',',num2str(current_point(2)),')']);
disp(['当前函数值：',num2str(f(current_point(1),current_point(2)))])
% 迭代更新
while true
    g = grad_f(current_point(1),current_point(2)); 
    H = hessian(current_point(1),current_point(2)); 
    % 计算搜索方向
    d = inv(H)*g;
    next_point = current_point - d;
    iter = iter + 1;
    disp(' ')
    disp(['迭代次数',num2str(iter),':']);
    disp(['当前点：（',num2str(next_point(1)),',',num2str(next_point(2)),')']);
    disp(['当前函数值：',num2str(f(next_point(1),next_point(2)))]);
    if norm(next_point - current_point) < precision
        break
    end
    current_point = next_point;
end
disp(' ')
disp(['迭代完成']);
disp(['最小值对应的点：（',num2str(next_point(1)),',',num2str(next_point(2)),')']);
disp(['函数最小值：',num2str(f(next_point(1),next_point(2)))]);
```

## 作业三-可视化无约束优化算法

### Q1

绘制$f(x,y)=x^4+2*y^4在x\in [-4,4],y\in [-4,4]$上的等高线

```matlab
% 绘制$f(x,y)=x^4+2*y^4在x\in [-4,4],y\in [-4,4]$上的等高线
f = @(x,y) x.^4 + 2*y.^4;
% 生成网格数据
[X,Y] = meshgrid(-4:0.1:4);
Z = f(X,Y)
% 绘制等高线图
figure
contour(X,Y,Z,'ShowText','on');
xlabel('x');
ylabel('y');
title('Contour Plot of f(x,y)=x^4+2y^4')
colorbar
```

<img src="C:\Users\keyang\AppData\Roaming\Typora\typora-user-images\image-20230531003934765.png" alt="image-20230531003934765" style="zoom:80%;" />



