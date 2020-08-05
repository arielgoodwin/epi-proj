% projecting onto the K-epigraph of F(X) = (1/2)XX^t
% K = S^n_+
% algorithm: subgradient method, based on algorithm 3.2 in Quang's writeup
% author: ariel goodwin

% Given: A pair (X_bar, V_bar) in the product space that we are trying to
% project onto the K-epigraph of F

% X_bar is a rectangular matrix (could be square)
% V_bar is a symmetric matrix
n = 2
m = 2

% pair to be projected
X_bar = rand(n,m);
V_bar = rand(n);
V_bar = V_bar*V_bar';

% initial point for the algorithm
V_o = rand(n);
V_o = V_o*V_o';
V_k = V_o
k = 1;

% iterate until termination condition
while k < 3000
    
    % step size parameter 
    L = max(1, norm((inv(V_k+eye(n))*X_bar),"fro"));
    
    % intermediate variables
    Y = V_k + V_bar;
    W = inv(V_k+eye(n));
    U = X_bar*X_bar';
    A = (1/(k*L))*(Y - 0.5*W*U*W);
    
    % update point and increment counter
    V_k = V_k - A;
    k = k+1;
    
end

% compute the projection 
proj1 = inv(V_k + eye(n))*X_bar
proj2 = V_k + V_bar