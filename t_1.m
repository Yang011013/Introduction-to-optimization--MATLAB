syms t;
f = t^4 - t^2 - 2*t + 5;
[x_optimization,f_optimization] = Golden_Selection_Method(f,-10,10);
disp([num2str(x_optimization),num2str(f_optimization)])

x = -2:0.001:2.5;
fx = x.^4 - x.^2 - 2.*x + 5;
figure(1)
plot(x,fx);