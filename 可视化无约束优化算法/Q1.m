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
