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