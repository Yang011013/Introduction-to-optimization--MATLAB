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