% 传参为3个长度为3的列向量，判断三个列向量是否线性相关，如果是则返回1，否则返回0
function flag = f5(a1,a2,a3)
    A = [a1; a2; a3]
    if det(A) == 0
        flag = 1
    else
        flag = 0
end