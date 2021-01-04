function theta = matrices2theta(Abar, Bbar, C, D, K, x0) 
    %%
    % Function INPUT
    % Abar      System matrix Abar (matrix of size n x n)
    % Bbar      System matrix Bbar (matrix of size n x m)
    % C         System matrix C (matrix of size l x n)
    % D         System matrix D (matrix of size l x m)
    % K         System matrix K (matrix of size n x l)
    % x0        Initial state (vector of size n x one)

    % Function OUTPUT
    % theta     Parameter vector (vector of size n*n+n*m+l*n+l*m+n*l+n)
    
    % Using full parameterization
    n = size(Abar, 1); m = size(Bbar, 2); l = size(C, 1);
    theta = zeros(n*n + n*m + l*n + l*m + n*l + n, 1);
    
    % A-matrix
    theta(1:n*n) = Abar;
    
    % B-matrix
    theta(n*n + 1 : n*n + n*m) = Bbar(:);
    
    % C-matrix
    theta(n*n + n*m + 1 : n*n + n*m + l*n) = C(:);
    
    % D-matrix
    theta(n*n + n*m + l*n + 1 : n*n + n*m + l*n + l*m) = D(:);

    % K-matrix
    theta(n*n + n*m + l*n + l*m + 1 : n*n + n*m + l*n + l*m + n*l) = K(:);

    % x0
    theta( n*n + n*m + l*n + l*m + n*l + 1 : end) = x0;
    
end