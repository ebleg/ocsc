function [A, B, C, D, x0, sv] = subspaceID(u, y, s, n, method)
    % Function INPUT 
    % u         system input (matrix of size N x m)
    % y         system output (matrix of size N x l)
    % s         block size (scalar)
    % n         model order (scalar)
    % method    method (string e.g. 'moesp')
    %
    % Function OUTPUT
    % A         System matrix A (matrix of size n x n)
    % B         System matrix B (matrix of size n x m)
    % C         System matrix C (matrix of size l x n)
    % D         System matrix D (matrix of size l x m)
    % x0        Initial state (vector of size n x one)
    % sv        Singular values (vector of size n x one)

    [N, l] = size(y);
    [~, m] = size(u);
    
    % Construct the hankel matricexs
    Y = hankel(y(1:l*s), y(l*s:end));
    U = hankel(u(1:m*s), u(m*s:end));
    
    Y2 = hankel(y(1:l*2*s), y(l*2*s:end));
    U2 = hankel(u(1:l*2*s), u(l*2*s:end));
    
    Up = U2(1:s,:);
    Uf = U2(s+1:end,:);
    Yp = Y2(1:s,:);
    Yf = Y2(s+1:end,:);   

    switch method
        case 'moesp'
            % LQ factorization 
            [~, L] = qr([U; Y]', 0); L = L';
            L22 = L(s*l+1:end, s*m+1:end); 

            % Singular value decomposition
            [U_svd, S_svd, ~] = svd(L22);
               
        case 'pi-moesp'
            % LQ factorization
            [~, L] = qr([Uf; Up; Yf]', 0); L = L';
            L32 = L(2*s+1:end, s+1:2*s); 

            % Singular value decomposition 
            [U_svd, S_svd, ~] = svd(L32);

        case 'po-moesp'
            % LQ factorization
            Z = [Up; Yp];
            [~, L] = qr([Uf; Z; Yf]', 0); L = L';
            L32 = L(3*s+1:end,s+1:2*s);

             % Singular value decomposition 
            [U_svd, S_svd, ~] = svd(L32);
    end
    
    % Output the sv's as vector instead of diag matrix
    sv = diag(S_svd);
    
    % Extract A and C from the extended observability matrix
    C   = U_svd(1:l, 1:n);
    A   = U_svd(1:s*l-l, 1:n)\U_svd(l+1:end, 1:n); 
    
    % Precompute & store CA^k for efficiency
    CA_powers = nan([size(C), N]);
    CA_powers(:,:,1) = C;
    
    for k = 2:(N+1)
        CA_powers(:,:,k) = CA_powers(:,:,k-1)*A;
    end

    % Preallocate the phi matrix for efficiency
    phi = zeros(l*N, n + n*m + l*m);

    % Construct the phi matrix
    for k = 0:size(y)-1
        % Build the convolution term
        conv_term = zeros(l,n); 
        for j = 0:k-1
              conv_term = conv_term + u(j+1)*CA_powers(:,:,k-j);
        end
        % Insert into the phi matrix
        phi(k*l+1:(k+1)*l, :) = [CA_powers(:,:,k+1) conv_term u(k+1)']; 
    end 
    
    % Extract B, x0 and D using least-squares
    tmp = pinv(phi)*y;
    x0 = tmp(1:n);
    B = reshape(tmp((n+1):(n + m*n)), [n, m]);
    D = reshape(tmp((n+m*n+1):end), [l, m]);
end