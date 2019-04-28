function X=HHL(f,x0,x1,tol,max)

This function uses the Aitken iterative algorithm to find all approximate zeros of the unary function f in a given interval

***************input value****************

f is the original function handle 

x0 is the left endpoint of the root interval (default is -10)

x1 is the right endpoint of the root interval (default is 10)

tol is the target error (default is 10^(-5))

max is the maximum number of iterations (default is 50)

%***************output value****************

X is the approximate solution of all zeros of the function in this interval

%**************Use example***************

1.>>myfun=@(x)x^3-6*x^2+11*x-6;

  >>X=HHL(myfun,0,5)

2. in the myfun.m file

   function y=myfun(x)

   y=x^3-6*x^2+11*x-6;

   >>X=HHL(@myfun,0,5)
