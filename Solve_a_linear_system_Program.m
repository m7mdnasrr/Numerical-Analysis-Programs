% Mohamed Nasr El-Din Mohamed Magdy Mohamed
% Numerical Analysis
% Solve a linear system Program

function Solve_a_linear_system_Program()
    disp('Solve a linear system AX = B');

    max_iter = 2004;

    while true
        A = input('Enter the coefficient matrix A: ');
        
        % Check if the input is valid
        if ismatrix(A) && size(A, 1) == size(A, 2)
            break;
        else
            disp('Invalid input. The coefficient matrix A must be a square matrix.');
        end
    end

    while true
        b = input('Enter the right-hand side vector B: ');
        
        % Check if the input is valid
        if isvector(b) && length(b) == size(A, 1)
            break;
        else
            disp('Invalid input. The vector B must have the same number of rows as matrix A.');
        end
    end

    % Automatically initialize the initial guess vector x0
    x0 = zeros(length(b), 1);

    while true
        tol = input('Enter the desired tolerance: ');
        
        % Check if the input is valid
        if isnumeric(tol) && isscalar(tol) && tol > 0
            break;
        else
            disp('Invalid input. The tolerance must be a positive number.');
        end
    end

    disp('Numerical Methods Available:');
    disp('(1) Jacobi Method');
    disp('(2) Gauss-Seidel Method');
    disp('(3) Relaxation Method');

    
    while true
        method = input('Choose the number of the method you want to use: ');
        
        % Check if the input is valid
        if isnumeric(method) && ismember(method, 1:3) && isscalar(method) && method == floor(method)
            break;
        else
            disp('Invalid method selected. Please try again.');
        end
    end

    if method == 3
        while true
            omega = input('Enter the relaxation factor (1 < omega < 2): ');
            
            % Check if the input is valid
            if isnumeric(omega) && isscalar(omega) && omega >= 1 && omega <= 2
                break;
            else
                disp('Invalid input. The relaxation factor omega must be between 1 and 2.');
            end
        end
    end

    switch method
        case 1 % Jacobi Method
            disp('Using Jacobi Method');
            [x, iter] = jacobi_method(A, b, x0, tol, max_iter);
        case 2 % Gauss-Seidel Method
            disp('Using Gauss-Seidel Method');
            [x, iter] = gauss_seidel_method(A, b, x0, tol, max_iter);
        case 3 % Relaxation Method
            disp('Using Relaxation Method');
            [x, iter] = relaxation_method(A, b, x0, tol, max_iter, omega);
    end

    fprintf('\nSolution reached in %d iterations:\n', iter);
    disp(x);
end

% Jacobi Method
function [x, iter] = jacobi_method(A, b, x0, tol, max_iter)
    n = length(b);
    x = x0;
    iter = 0;
     fprintf(' %-6s', 'k');
    for i = 1:n
        fprintf(' %-11s', ['x', num2str(i), '^(k)']);
    end
    fprintf('%-10s\n', '|x^(k) - x1^(k-1)|');
    fprintf('%s\n', repmat('-', 1, 8 + 12*n + 18)); % Separator line
    fprintf(' %-6d', iter); % Initial guess as iteration 0
    fprintf('%-12.6f', x);
    fprintf('  %-12.6f\n', norm(x, inf));
    while iter < max_iter
        x_new = zeros(n, 1);
        for i = 1:n
            x_new(i) = (b(i) - A(i, [1:i-1, i+1:n]) * x([1:i-1, i+1:n])) / A(i, i);
        end
        iter = iter + 1;
        acc = norm(x_new - x, inf);
        fprintf(' %-6d', iter);
        fprintf('%-12.6f', x_new);
        fprintf('  %-12.6f\n', acc);
        x = x_new;
        if acc < tol
            break;
        end
    end
end

% Gauss-Seidel Method
function [x, iter] = gauss_seidel_method(A, b, x0, tol, max_iter)
    n = length(b);
    x = x0;
    iter = 0;
    fprintf(' %-6s', 'k');
    for i = 1:n
        fprintf(' %-11s', ['x', num2str(i), '^(k)']);
    end
    fprintf('%-10s\n', '|x^(k) - x1^(k-1)|');
    fprintf('%s\n', repmat('-', 1, 8 + 12*n + 18)); % Separator line
    fprintf(' %-6d', iter); % Initial guess as iteration 0
    fprintf('%-12.6f', x);
    fprintf(' %-12.6f\n', norm(x, inf));
    while iter < max_iter
        x_new = x;
        for i = 1:n
            x_new(i) = (b(i) - A(i, [1:i-1, i+1:n]) * x_new([1:i-1, i+1:n])) / A(i, i);
        end
        iter = iter + 1;
        acc = norm(x_new - x, inf);
        fprintf(' %-6d', iter);
        fprintf('%-12.6f', x_new);
        fprintf(' %-12.6f\n', acc);
        x = x_new;
        if acc < tol
            break;
        end
    end
end

% Relaxation Method
function [x, iter] = relaxation_method(A, b, x0, tol, max_iter, omega)
    % Normalize matrix A to make diagonal elements 1
    [A, b] = normalize_diagonal(A, b);

    % Initialization
    n = length(b);          
    x = x0;                
    iter = 0;               
    r = b - A * x;          

    fprintf(' %-6s', 'k');
    for i = 1:n
        fprintf(' %-11s', ['x', num2str(i), '^(k)']);
    end
    for i = 1:n
        fprintf(' %-11s', ['r', num2str(i), '^(k)']);
    end
    fprintf('\n%s\n', repmat('-', 1, 10 + 23 * n)); % Separator line

    % Initial values (iteration 0)
    fprintf(' %-6d', iter);
    fprintf('%-12.6f', x);
    fprintf('%-12.6f', r);
    fprintf('\n');

    while iter < max_iter
        x_old = x;          
        
        % Find the index with the largest residual
        [~, idx] = max(abs(r));

        % Update only the variable with the largest residual
        sigma = 0;          
        for j = 1:n
            if j ~= idx
                sigma = sigma + A(idx, j) * x(j);
            end
        end
        
        x(idx) = (1 - omega) * x_old(idx) + omega * (b(idx) - sigma) / A(idx, idx);
        
        r = b - A * x;

        iter = iter + 1;

        fprintf(' %-6d', iter);
        fprintf('%-12.6f', x);
        fprintf('%-12.6f', r);
        fprintf('\n');

        % Check convergence
        if norm(r, inf) <= tol
            break;
        end
    end
end

% Function to normalize the diagonal of A to 1
function [A, b] = normalize_diagonal(A, b)
    n = size(A, 1); % Number of rows
    for i = 1:n
        factor = A(i, i);  % Factor to make the diagonal 1
        A(i, :) = A(i, :) / factor;
        b(i) = b(i) / factor;
    end
end
