% 传参为一个3*3的矩阵，返回该矩阵的逆矩阵，如果无逆矩阵，返回-1
function flag = f3(A)
    if det(A) == 0 % 行列式等于0，矩阵不可逆
        flag = -1
    else
        flag = inv(A)
    end
end