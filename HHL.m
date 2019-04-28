function X=HHL(f,x0,x1,tol,max)

% This function uses the Aitken iterative algorithm to find all approximate zeros of the unary function f in a given interval
%X=HHL(f,x0,x1,tol,max)
%***************input value****************
%f is the original function handle 
%x0 is the left endpoint of the root interval (default is -10)
%x1 is the right endpoint of the root interval (default is 10)
%tol is the target error (default is 10^(-5))
%max is the maximum number of iterations (default is 50)
%***************output value****************
%X is the approximate solution of all zeros of the function in this interval
%**************Use example***************
%1.>>myfun=@(x)x^3-6*x^2+11*x-6;
%  >>X=HHL(myfun,0,5)
%2. in the myfun.m file
%   function y=myfun(x)
%   y=x^3-6*x^2+11*x-6;
%   >>X=HHL(@myfun,0,5)

%**************Local variable**************
%i loop variable
%j loop variable
%k current function zeros +1
%xk iteration initial value
%step iteration initial value step
%flag repeat flag variable, =0 no repeat, =1 repeat
%Q  ~x in the formula
%R  ￣x n the formula
%P value after each iteration
%x The final iteration result of the initial value of each iteration
%swap_temp data exchange staging variable
%str string temporary variable
%temp structure temporary variables

%**************************************Input variable processing*****************************************
if(nargin<5)
    max=50;   % is not equal to 50 when the maximum number of iterations is entered
end
if(nargin<4)
    tol=10^(-5);   % is equal to 10^(-5) when there is no input error
end
if(nargin==2) 
    error('[Please enter the third parameter] can not only enter the left or right endpoint, you must enter or use the default interval at the same time. For details, please refer to help HHL'); % can not only enter the left endpoint, not the right endpoint
end
if(nargin==1)
    x0=-10;        % is not equal to [-10,10] when inputting the root interval
    x1=10;
end
if(nargin<1)
    error('[Please enter the first parameter, that is, the function handle] must input the original function, please re-enter. For details, please refer to help HHL'); % must input the original function handle, no original function can not be solved
end

%对输入的区间查错
if(x1<=x0)
    error('[The third parameter should not be smaller than the second parameter] The right end of the interval is less than or equal to the left end point, and the interval is incorrect. Please re-enter. For details, please refer to help HHL');
end
if(x1-x0>200000)    % prevents the solution interval from being too 
    error('[The third parameter should not be much larger than the second parameter] The length of the interval is too large, and the solution requires a lot of time. Please reduce the length of the interval. For details, please refer to help HHL');
end
if(max>1000||max<0)   
    error('[The fifth parameter should not be too large or negative] The number of iterations is unreasonable, please re-enter. For details, please refer to help HHL');
end
if(tol<0)
    error('[The fourth parameter should not be negative] The error limit should be positive, please re-enter. For details, please refer to help HHL');
end
%**************************************initialization***************************************
k=1;  %current function has zero number of zeros
step=(x1-x0)/1000;  %iteration initial value step
if(step>0.3)        %prevents the iteration initial value step size from being too large, resulting in a leak root
    step=0.3;
end

%Determine if the special point is root
if(abs(feval(f,x0))<tol)  %Check if the left endpoint is root
    X(k)=x0;
    k=k+1;
end
if(abs(feval(f,x1))<tol)  %Check if the right endpoint is root
    X(k)=x1;
    k=k+1;
end
if(x0<0&&x1>0&&abs(feval(f,0))<tol)   %Check if the right endpoint is root
    X(k)=0;
    k=k+1;
end

%************************Function core, repeatedly using Aitken algorithm to calculate the result******************************
for xk=x0:step:x1;  %takes a large number of iteration initial values
    % first iteration
    Q=feval(f,xk)+xk;   % Substituting the initial value into the iterative formula. Iterative formula: f(x)+x, the same below
    R=feval(f,Q)+Q;
    if((R-2*Q+xk)==0)    %P(1) The denominator is 0, the iterative result does not converge, and the iteration of the initial value of the next iteration is directly performed.
        continue;
    end    
    P(1)=R-((R-Q)^2)/(R-2*Q+xk);  % first iteration result
    
    % loop iteration until the difference between two iterations is less than the error limit
    for j=2:max;
        Q=feval(f,P(j-1))+P(j-1);
        R=feval(f,Q)+Q;
        if((R-2*Q+P(j-1))==0)  %P(j) denominator is 0, the iterative result does not converge, and the iteration of the initial value of the next iteration is directly performed.
            break;
        end
        P(j)=R-((R-Q)^2)/(R-2*Q+P(j-1));  % jth iteration result
        err=abs(P(j)-P(j-1));  % current error
        x=P(j);
        if (err<tol) % The absolute value of the difference between two iterations is less than the end of the error limit
            break;
        end
    end
    
%******************************************Judge whether the value iterated out is Compliance with the conditions of the root **************************************************
    flag=0;   %repeat flag variable initial value
    if(exist('x','var')==0)   %does not find x, does not have to judge the root, directly iterate the initial value of the next iteration
        continue;
    end
    if(j<max &&x>x0&&x<x1 &&abs(feval(f,x))<tol)  %Root requires conditions: 1. Iteration results converge. 2. The iterative result is within the root interval. 3. The iteration result is the zero point of the function
        for i=1:k-1   
            if(abs(X(i)-x)<10^(-3))  %Check whether the result of this iteration is repeated with the existing root
                flag=1;
                break;
            end
        end
        if(flag==0)   %is not repeated
            X(k)=x;   %is not repeated
            k=k+1;
        end
    end
end

%***********Bulk sorting of roots, from small to large array ******************
for i=1:k-1;
    for j=i:k-1;
        if(X(i)>X(j))   %exchanges two data
            swap_temp = X(j);
            X(j) = X(i);
            X(i) = swap_temp;
        end
    end
end

%**************************Zero situation report ****************************
%get the original function
temp=functions(f);  %Get function handle details
str=temp.function;  % get the function handle content
if(strcmp(temp.type,'anonymous'))  %is an anonymous function?    
    [i,j]=size(str); 
    for i=1:j;
        if(str(i)==')')  %search for the character ‘)’
            break;
        end
    end
    str=str(i+1:end);  %takes the characters after ‘)’. For example, "fun=@(x)x^2-1" is intercepted as "x^2-1"
    
else   % is the m file function
    str=strcat(str,'.m');
    fp = fopen(str,'r');  %Open the corresponding m file
    str=fscanf(fp,'%s');  %Read m file content
    [i,j]=size(str);   %Get the length of the string
    flag=0;  % flag variable initial value
    for i=1:j;
        if(str(i)=='=')  %Get the length of the string
            if(flag==1)  % second search for ‘=’
                 break;
            end   
            flag=1;  % second search for ‘=’
        end
    end
    str=str(i+1:end-1);  % takes the second ‘=’ character. For example, "function y=fun(x)y=x^2-1;" is intercepted as "x^2-1"  
end
str=strcat(str,'=0'); % adds "=0" to the string obtained above

if(k==1)  % adds "=0" to the string obtained above
    ('equation %s has no approximate solution \n' in the interval [%g, %g], str, x0, x1);
else   has root condition
    fprintf('equation %s has a total of %d approximate solutions in the interval [%g, %g], arranged from small to large as follows: ', str, x0, x1, k-1);
end
