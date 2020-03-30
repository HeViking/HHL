# 函数简介
`function X = HHL(f,x0,x1,tol,max) `  
本函数用于求出一元若干次函数f在所给定区间中的所有近似零点  
用Matlab2014实现

# 具体介绍
  
#### 输入值
`f` 为函数句柄   
`x0` 为求根区间左端点(默认值为-10)  
`x1` 为求根区间右端点(默认值为10)  
`tol` 为目标误差(默认值为 10^(-5) )  
`max` 为最大迭代次数(默认值为50)  
  
#### 输出值
`X` 为该函数在该区间里所有零点的近似解  
  
#### 使用范例
1.在命令行中输入函数  
`myfun = @(x)x^3 - 6*x^2 + 11*x - 6;`  
`X = HHL(myfun,0,5)`  

或者
2.在myfun.m文件中输入函数  
`function y=myfun(x)`  
`y = x^3 - 6*x^2 + 11*x - 6;`  
在命令行中     
`X = HHL(@myfun,0,5)`  

输出结果如下：  
方程 x^3-6*x^2+11*x-6=0 在区间\[0,5\]有3个近似解，从小到大排列如下：  
X = 1.0000    2.0000    3.0000
