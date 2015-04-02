function [x] = gauss_seidel(A, b, x0, iters)
    n = length(A);
    x = x0;
    for k = 1:iters
        for i = 1:n
            x(i) = (1/A(i, i))*(b(i) - A(i, 1:n)*x + A(i, i)*x(i));
        end
    end
end
