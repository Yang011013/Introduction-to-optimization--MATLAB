% 传参为一个5*5的矩阵，判断该矩阵是否为对称矩阵，如果是则返回1，否则返回0
function flag = f4(A)
    if A' == A
        flag = 1
    else
        flag = 0      
end