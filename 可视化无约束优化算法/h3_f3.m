% 定义函数和梯度
f = @(x,y) x.^4 + 2*y.^4;
grad_f = @(x,y) [4*x^3; 8*y^3];

% 初始化
current_point = [4; 4];
precision = 0.1;

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
hold on

% 迭代更新
while true
    gradient = grad_f(current_point(1),current_point(2));
    % 绘制迭代点
    plot(current_point(1),current_point(2),'bo','MarkerSize',4,'LineWidth',1);
    text(current_point(1)-0.2,current_point(2)-0.1,['(',num2str(current_point(1)),',',num2str(current_point(2)),')']);
    text(current_point(1)-0.2,current_point(2)-0.3,['f=',num2str(f(current_point(1),current_point(2)))]);
    
    % 精确搜索, 找出[-100,100]范围内的alpha取值
    alpha = fminbnd(@(alpha)f(current_point(1)-alpha*gradient(1),current_point(2)-alpha*gradient(2)),-100,100);
    next_point = current_point - alpha * gradient;
    
    % 绘制迭代箭头
    arrow_vec = next_point - current_point;
    quiver(current_point(1),current_point(2),arrow_vec(1),arrow_vec(2),'Color','green','LineWidth',1);
    
    if norm(next_point - current_point) < precision
        break
    end
    current_point = next_point;    
end
% 绘制最小点
plot(current_point(1),current_point(2),'ro','MarkerSize',4,'LineWidth',2);
text(current_point(1),current_point(2)+0.2,['(',num2str(current_point(1)),',',num2str(current_point(2)),')']);
text(current_point(1),current_point(2)+0.4,['f=',num2str(f(current_point(1),current_point(2)))]);
disp('Optimization complete.');
disp(['Minimum point:(', num2str(current_point(1)),',',num2str(current_point(2)),')'])
disp(['Minimum value:',num2str(f(current_point(1),current_point(2)))]);
