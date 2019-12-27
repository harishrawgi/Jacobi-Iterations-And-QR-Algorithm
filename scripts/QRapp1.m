%% Function description
% The function calculates eigen vectors and eigen values of a matrix using
% an iterative method
% Inputs: A (a matrix)
%         N (max iterations)
%         ch (choice for gs or mgs), 0 means gs and 1 means mgs
% Outputs: evalues (eigen values for A)
%          evectors (eigen vectors for A)
function [evalues, evectors] = QRapp1(A, N, ch)
    
    % get the size of the matrix
    n = size(A,1);
    
    % initialize the eigen vector matrix
    evectors = eye(n);
    
    % set A as the first iterate
    A_k = A;
    
    % iteratively use QR decomposition to compute eigen values and vectors
    for i = 1:N
        
        % use gs or mgs depending on the choice supplied
        if (ch == 0)
            [Q_k, R_k] = gs(A_k);
        else
            [Q_k, R_k] = mgs(A_k);
        end
        
        % update eigen vectors for the current iterations
        evectors = evectors*Q_k;
        
        % compute the next iterate using R and Q
        A_k = R_k*Q_k;
    end
    
    % the eigen values are the diagonal entries of the last iterate
    evalues = diag(A_k);
end