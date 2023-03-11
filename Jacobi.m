%% Notes
% please note that square matrix of coefficients must  
% be diagonally dominant for this program to work.
%% Initializing
%A = input('Please enter the coefficient matrix:\n') ;
%b = input('Please enter the source vector:\n') ;
A = [4 2 1;1 6 1;0 3 -5];
b = [3;7;-2];
n = length(b); % number of equations
x = zeros(n,1); % the solution vector
xnew = zeros(n,1); % the xi* vector
% by Jacobi's method algorithm:
% xi*(in each iteration) = 
% [-1/a(i,i)]*[sum(ai,j*xj) - bi] for j = 1:n & j~=i
% if j == i then we will have repeated values on rhs and lhs
x(:) = 0; % the guess/starter vector
iterlimit = 10; % the maximum number of iterations
tol = 0.001; % the absolute value of the difference between x and xnew 
% if abs(x-xnew)>tol then we need more iterations
%% Jacobi's method into codes
for iteration = 1:iterlimit % loop of iterations
    convergence = true; % boolean value
    for i=1:n % the loop of equations
        Sum = 0; % the values of previous equations will not be added to the current summation
        for j = 1:n % the loop of summation
            if j~=i
                Sum = Sum + A(i,j)*x(j); % according to the algorithm
            end
        end
        xnew(i) = -1/A(i,i) * (Sum - b(i)); % according to the algorithm
        if abs(xnew(i) - x(i)) > tol % we don't have the convergance
            convergence = false;
        end
    end
    if convergence
        break ; % break the loop of iterations
    end
    x = xnew;
    
    disp('iterations:') ;
    iteration
    disp('solution:')
    xnew
end