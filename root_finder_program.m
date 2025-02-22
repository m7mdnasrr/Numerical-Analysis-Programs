% God Bless Mohamed Nasr El-Din Mohamed Magdy Mohamed
% Numerical Analysis
% Root Finder Program

function root_finder_program()
    disp('Root Finder Program:');
    
    maxIterations = 2004;
    
    func = input('Enter the function F(x)= ', 's');
    f = str2func(['@(x) ' func]);
    
    disp('Numerical Methods Available:');
    disp('(1) Bisection Method');
    disp('(2) False Position Method');
    disp('(3) Simple Fixed-Point Iteration Method');
    disp('(4) Newton-Raphson Method');
    disp('(5) Secant Method');
    
    while true
        try
            method = input('Choose the number of the method you want to use: ');
            
            % Check if the input is valid
            if isnumeric(method) && ismember(method, 1:5) && isscalar(method) && method == floor(method) 
                break;
            else
                disp('Invalid method selected, Please try again.');
            end
    end
end

    switch method
        case 1 % Bisection Method
            a = input('Enter the lower boundry of the interval: ');
            b = input('Enter the upper boundry of the interval: ');
            if f(a) * f(b) > 0
                disp('Error: Function must have opposite signs at a and b.');
                return;
            end
            tolerance = input('Enter the tolerance: ');
            [root, found] = bisectionMethod(f, a, b, tolerance, maxIterations);
            
        case 2 % False Position Method
            a = input('Enter the lower boundry of the interval: ');
            b = input('Enter the upper boundry of the interval: ');
            if f(a) * f(b) > 0
                disp('Error: Function must have opposite signs at a and b.');
                return;
            end
            tolerance = input('Enter the tolerance: ');
            [root, found] = falsePositionMethod(f, a, b, tolerance, maxIterations);
            
        case 3 % Simple Fixed-Point Iteration Method
            g_str = input('Enter the simple fixed-point function g(x): ', 's');
            g = str2func(['@(x) ' g_str]);
            initialGuess = input('Enter the initial guess: ');
            tolerance = input('Enter the tolerance: ');
            [root, found] = fixedPointMethod(g, initialGuess, tolerance, maxIterations);
            
        case 4 % Newton-Raphson Method
            df_str = input('Enter the derivative of the function: ', 's');
            df = str2func(['@(x) ' df_str]);
            initialGuess = input('Enter the initial guess: ');
            tolerance = input('Enter the tolerance: ');
            [root, found] = newtonRaphsonMethod(f, df, initialGuess, tolerance, maxIterations);
            
        case 5 % Secant Method
            initialGuess1 = input('Enter the first initial guess: ');
            initialGuess2 = input('Enter the second initial guess: ');
            tolerance = input('Enter the tolerance: ');
            [root, found] = secantMethod(f, initialGuess1, initialGuess2, tolerance, maxIterations);
            
        otherwise
            disp('Invalid method selected');
            return;
    end
    
    if found
        fprintf('Root found = %f\n', root);
    else
        disp('Method did not converge within the maximum number of iterations.');
    end
end

% Bisection Method
function [root, found] = bisectionMethod(f, a, b, tol, maxIterations)
    found = false;
    fprintf(' %-5s%-10s%-10s%-10s%-10s\n', ' i', '    a', '    b', '    c', '    f(c)');
    fprintf('%s\n', repmat('-', 47, 1)); % Line separator
    for iter = 0:maxIterations
        c = (a + b) / 2;
        fprintf('  %-5d%-10.6f%-10.6f%-10.6f%-10.6f\n', iter, a, b, c, f(c));
        if abs(f(c)) < tol || abs(b - a) < 2*tol
            root = c;
            found = true;
            return;
        elseif f(a) * f(c) < 0
            b = c;
        else
            a = c;
        end
    end
    root = (a + b) / 2;
end

% False Position Method
function [root, found] = falsePositionMethod(f, a, b, tol, maxIterations)
    found = false;    
    fprintf(' %-5s%-10s%-10s%-10s%-10s\n', ' i', '    a', '    b', '    c', '    f(c)');
    fprintf('%s\n', repmat('-', 47, 1)); % Line separator
    for iter = 0:maxIterations
        c = (a * f(b) - b * f(a)) / (f(b) - f(a));
        fprintf('  %-5d%-10.6f%-10.6f%-10.6f%-10.6f\n', iter, a, b, c, f(c));
        if abs(f(c)) < tol || abs(b - a) < tol
            root = c;
            found = true;
            return;
        elseif f(a) * f(c) < 0
            b = c;
        else
            a = c;
        end
    end
    root = (a * f(b) - b * f(a)) / (f(b) - f(a));
end

% Simple Fixed-Point Iteration Method
function [root, found] = fixedPointMethod(g, x0, tol, maxIterations)
    found = false;    
    fprintf(' %-5s%-10s%-10s%-10s\n', ' i', '   x0', '   x1',  ' |x1 - x0|');
    fprintf('%s\n', repmat('-', 36, 1)); % Line separator
    for iter = 0:maxIterations
        x1 = g(x0);
        fprintf('  %-5d%-10.6f%-10.6f%-10.6f\n', iter, x0, x1, abs(x1 - x0));
        if abs(x1 - x0) <= tol
            root = x1;
            found = true;
            return;
        end
        x0 = x1;
    end
    root = x0;
end

% Newton-Raphson Method
function [root, found] = newtonRaphsonMethod(f, df, x0, tol, maxIterations)
    found = false;
    fprintf(' %-5s%-10s%-10s%-10s\n', ' i', '   x0', '   x1',  ' |x1 - x0|');
    fprintf('%s\n', repmat('-', 36, 1)); % Line separator
    for iter = 0:maxIterations
        if abs(df(x0)) < eps
            disp('Derivative near zero, Newton-Raphson may fail.');
            root = NaN;
            return;
        end
        x1 = x0 - f(x0) / df(x0);
        fprintf('  %-5d%-10.6f%-10.6f%-10.6f\n', iter, x0, x1, abs(x1 - x0));
        if abs(x1 - x0) <= tol
            root = x1;
            found = true;
            return;
        end
        x0 = x1;
    end
    root = x0;
end

% Secant Method
function [root, found] = secantMethod(f, x0, x1, tol, maxIterations)
    found = false;
    fprintf(' %-5s%-10s%-10s%-10s%-10s\n', ' i', '   x0', '    x1', '   x2', ' |x2 - x1|');
    fprintf('%s\n', repmat('-', 47, 1)); % Line separator
    for iter = 0:maxIterations
        if abs(f(x1) - f(x0)) < eps
            disp('Division by near-zero value, Secant method may fail.');
            root = NaN;
            return;
        end
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        fprintf('  %-5d%-10.6f%-10.6f%-10.6f%-10.6f\n', iter, x0, x1, x2, abs(x2 - x1));
        if abs(x2 - x1) <= tol
            root = x2;
            found = true;
            return;
        end
        x0 = x1;
        x1 = x2;
    end
    root = x1;
end
