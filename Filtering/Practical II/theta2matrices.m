function [Abar, Bbar, C, D, K, x0] = theta2matrices(theta, Asize, Bsize, ...
                                                Csize, Dsize, Ksize, xsize)
    %%
    % Function INPUT
    % theta     Paramter vector (vector of size n*n+n*m+l*n+l*m+n*l+n)
    % Asize     Size of Abar 
    % Bsize     Size of Bbar 
    % Csize     Size of C 
    % Dsize     Size of D 
    % Ksize     Size of K 
    % xsize     Size of x0

    % Function OUTPUT
    % Abar      System matrix A (matrix of size n x n)
    % Bbar      System matrix B (matrix of size n x m)
    % C         System matrix C (matrix of size l x n)
    % D         System matrix D (matrix of size l x m)
    % K         System matrix K (matrix of size n x l)
    % x0        Initial state (vector of size n x one)
    
    % Use the output normal form parametrization
    n = Asize(1);
    m = Bsize(2);
    l = Csize(1);
    
    % Caution: SISO only!
    % A-matrix
    Abar = reshape(theta(1:n*n), Asize);
    
    % B-matrix
    Bbar = reshape(theta(n*n + 1 : n*n + n*m), Bsize);
    
    % C-matrix
    C = reshape(theta(n*n + n*m + 1 : n*n + n*m + l*n), Csize);
    
    % D-matrix
    D = reshape(theta(n*n + n*m + l*n + 1 : n*n + n*m + l*n + l*m), Dsize);

    % K-matrix
    K = reshape(theta(n*n + n*m + l*n + l*m + 1 : n*n + n*m + l*n + l*m + n*l), Ksize);

    % x0
    x0 = theta( n*n + n*m + l*n + l*m + n*l + 1 : end);
end