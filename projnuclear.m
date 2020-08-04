% Projection onto the nuclear norm unit ball
% author: ariel goodwin

% Y is a matrix, not necessarily square, stopThr is the termination
% condition for convergence, and guess is the initial guess lambda_o for
% the semismooth Newton method
function proj = projnuclear(Y,stopThr,guess)

    % perform singular value decomposition
    [U,S,V] = svd(Y);

    % check if the vector lies in the ball already
    if sum(diag(S,0)) <= 1
        proj=Y;
        return
    end
    
    % termination condition
    delta = stopThr;
    
    % initial guess and iteration counter
    lambda_k = guess;
    k = 0;

    % isolate the singular values
    sing = diag(S,0);

    % Iterate while checking termination condition
    
    while abs(F(lambda_k,sing)) >= delta
    
        % determine an element of the B-subdifferential
        g_k = sum(sing - lambda_k,0);

        % find the descent direction (if g_k = 0, d_k = -inf so d_k
        % will then become -lambda_k)
        d_k = -F(lambda_k,sing)/(g_k);
        d_k = max(-lambda_k,d_k);

        % update test point and increment index
        lambda_k = lambda_k + d_k;
        k=k+1;
    
    end

    % compute projection
    d = max(sing - lambda_k,0);
    D = diag(d);
    proj = U*D*V';
  

end
%% 
% Auxilliary function: grad = derivative of Theta w.r.t. lambda where Theta 
% is the function from Hoheisel's paper

function grad = F(lambda,sing)
   grad = 1 - sum(max(sing-lambda,0));
end