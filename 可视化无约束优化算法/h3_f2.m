% 绘制$f(x,y)=x^4+2*y^4在x\in [-4,4],y\in [-4,4]$上的三维图
f = @(x,y) x.^4 + 2*y.^4;
% 生成网格数据
[X,Y] = meshgrid(-4:0.1:4);
Z = f(X,Y)
% 绘制三维图
figure
surf(X,Y,Z);
xlabel('x');
ylabel('y');
title('3D Plot of f(x,y)=x^4+2y^4')