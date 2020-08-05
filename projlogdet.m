% Projection onto the level sets of negative logdet function
% author: ariel goodwin

% X is a square matrix, stopThr is the termination
% condition for convergence, and guess is the initial guess lambda_o for
% the semismooth Newton method
function proj = projlogdet(X,stopThr,guess)

    % check if the matrix lies in the level set
    % note that this corresponds to alpha = 0, the 0-level set for the
    % function. 
    % in general, these correspond to the condition detX >= exp(-alpha)
    % for the alpha-level set
    if det(X) >= 1 
        
        % it is in the level set
        proj = X;
        return;
        
    end
       
    % termination condition
    delta = stopThr;
    
    % initial guess and iteration counter
    lambda_k = guess;
    k = 0;
    
    % perform a singular value decomp
    [U,S,V] = svd(X);
    eg = eig(X);

    % Iterate while checking termination condition
    while abs(F(lambda_k,eg)) >= delta
    
        % now to compute the B-subdifferential of F at lambda_k
        g_k = G(lambda_k,eg);
        
        % small positive sequence going to zero
        epsilon_k = 1/2^(k+1);

        % find the descent direction
        d_k = -F(lambda_k,eg)/(g_k);
        d_k = max(-lambda_k + epsilon_k,d_k);
        
        % update test point and increment index
        lambda_k = lambda_k + d_k;
        k=k+1;
    
    end

    % compute projection
    d = 0.5*(eg+sqrt(eg.*eg + 4*lambda_k));
    D = diag(d);
    proj = U*D*V';

end
%% 
% %Function definitions: these are the first and second derivatives, respectively, 
% of Theta defined in Hoheisel's paper, w.r.t. lambda

function grad = F(lambda,v)

   u = 0.5*(v + sqrt(v.*v + 4*lambda));
   grad = sum(log(u));
end

function gradgrad = G(lambda,x)

    v = x.*x+4*lambda+x.*sqrt(x.*x+4*lambda);
    v = 1./v;
    gradgrad = 2*sum(v);
    
end