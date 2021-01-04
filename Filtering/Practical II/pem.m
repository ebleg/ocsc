function [Abar, Bbar, C, D, K, x0] = pem(A0, B0, C0, D0, K0, x00, y, ...
                                                               u, maxiter)
    %% Instructions:

    % Function INPUT 
    % A0        Initial guess for system matrix A (matrix of size n x n)
    % B0        Initial guess for system matrix B (matrix of size n x m)
    % C0        Initial guess for system matrix C (matrix of size l x n)
    % D0        Initial guess for system matrix D (matrix of size l x m)
    % K0        Initial guess for system matrix K (matrix of size n x l)
    % x00       Initial guess for initial state (vector of size n x one)
    % u         System input (matrix of size N x m)
    % y         System output (matrix of size N x l)
    % maxiter   Maximum number of iterations (scalar)
    %
    % Function OUTPUT
    % Abar      Estimate of system matrix A (matrix of size n x n)
    % Bbar      Estimate of system matrix B (matrix of size n x m)
    % C         Estimate of system matrix C (matrix of size l x n)
    % D         Estimate of system matrix D (matrix of size l x m)
    % K         Estimate of system matrix K (matrix of size n x l)
    % x0        Estimate of initial state (vector of size n x one)
    
    %% Prepare Gauss-Newton
    % Define matrix sizes
    Asize = size(A0); Bsize = size(B0); Csize = size(C0); Dsize = size(D0);
    xsize = size(x00); Ksize = size(K0);
    n = Asize(1);
    m = Bsize(2);
    l = Csize(1);
    
    % Initial theta
    theta = matrices2theta(A0, B0, C0, D0, K0, x00);
    p = length(theta);
    N = length(y);
    
    % Precompute the matrix derivatives with finite differences
    % (there are even more trivial ways to do this)
    h = 1e-6; % FD step size
    
    % Initialize matrices
    dAdth = zeros([Asize, p]);
    dBdth = zeros([Bsize, p]);
    dCdth = zeros([Csize, p]);
    dDdth = zeros([Dsize, p]);
    dKdth = zeros([Ksize, p]);
    dx0dth = zeros([xsize, p]);
    
    for i = 1:p
        % Dinite-difference vector
        h_vec = zeros(p, 1); h_vec(3) = h;
        [A_dist, B_dist, C_dist, D_dist, K_dist, x0_dist] = ...
                theta2matrices(theta + h_vec, Asize, Bsize, Csize, ...
                                                     Dsize, Ksize, xsize);
        % Finite difference computation
        dAdth(:,:,i) = (A_dist - A0)/h;                                         
        dBdth(:,:,i) = (B_dist - B0)/h;                                         
        dCdth(:,:,i) = (C_dist - C0)/h;                                         
        dDdth(:,:,i) = (D_dist - D0)/h;                                         
        dKdth(:,:,i) = (K_dist - K0)/h;                                         
        dx0dth(:,:,i) = (x0_dist - x00)/h;                                         
    end
    
    %% Gauss-Newton
    % Set convergence checks to false
    converged = false;
    maxiter_reached = false;
    step_size = 1;
    lambda = 0.5;
    iter = 1;
    
    while ~converged && ~maxiter_reached
        % Find the matrices based on current theta
        [A, B, C, D, K, x0] = theta2matrices(theta, Asize, Bsize, ...
                                              Csize, Dsize, Ksize, xsize);                                          
        
        % Compute E
        [y_sim, x_sim] = simsystem(A, [B K], C, [D zeros(l, l)], x0, [u y]);
        E = y - y_sim;
        
        % Compute phi
        phi = zeros(l*N, p);
        for i = 1:p
            % Simulate 'partial derivative' system
            dydth = simsystem(A, [dAdth(:,:,i) dBdth(:,:,i) dKdth(:,:,i)], ...
                C, [dCdth(:,:,i), dDdth(:,:,i), zeros(l, l)], ...
                dx0dth(:,:,i), [x_sim, u, y_sim]);            
            phi(:, i) = -dydth;
        end
        
        % Perform time step
        J = 2/N*phi'*E;
        H = 2/N*(phi'*phi);
        
        theta_new = theta - step_size*(H + eye(size(H))*lambda)\J;
        
        % Check convergence
        if norm(theta_new - theta) < 1e-8
            converged = true;
        elseif iter == maxiter
            maxiter_reached = true;
            warning('Maximum iterations reached');
        end
        theta = theta_new;
        iter = iter + 1;
    end
   [Abar, Bbar, C, D, K, x0] = theta2matrices(theta, Asize, Bsize, ...
                                               Csize, Dsize, Ksize, xsize);
end