%% Function description
% The function solves a system of equations for 
% diagonally dominant matrix A using the Jacobi iteration method
% Inputs: A (diagonally dominant matrix)
%         b (a vector)
%         tol (specifies the tolerance of accuracy)

%% Function code
function [x, iteration_count] = jacobi(A, b, tol)

% get the size of a
n = size(A,1);

% initialize the iteration count
iteration_count = 0;

% intitialize the beignning point
x_0 = zeros(n,1);

% calculate the initial residue
r_0 = A*x_0 - b;

% set x_k as x_0
x_k = x_0;

% initialize x_kplus1 for efficiency of code
x_kplus1 = zeros(n,1);

relative_residue = 1;

% iteratively calculate x using Jacobi method
while relative_residue > tol
    
    % perform one iteration of Jacobi method
    for i=1:n
        
        % initialize sum of multiplied terms to zero
        sigma = 0;
        
        % iterate over all terms and multiply accordingly and add them
        for j=1:n
            if j~=i
                sigma = sigma + A(i,j)*x_k(j);
            end
        end
        
        % update x as per the iteration step
        x_kplus1(i) = (b(i) - sigma)/A(i,i);
    end
    
    %update the iteration count
    iteration_count = iteration_count + 1;
    
    % compute the new residual
    r_kplus1 = A*x_kplus1 - b;
    
    % calculate the relative residue
    relative_residue = norm(r_kplus1)/norm(r_0);
    
    % Break the iteration if we reach the required tolerance
    %if(relative_residue < tol)
    %    break;
    %end
    
    % make the new x as the old x for the next iteration
    x_k = x_kplus1;
    
end

% set x to the last iterate
x = x_kplus1;

end