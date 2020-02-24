% % Non-negative Matrix Factorization via Nesterov's Optimal Gradient Method: Improved via Randomized Subspace Iterations NeNMF (RSI) 

% Reference
%  F. Yahaya, M. Puigt, G. Delmaire, G. Roussel, Faster-than-fast NMF using
%  random projections and Nesterov iterations, to appear in the Proceedings
%  of iTWIST: international Traveling Workshop on Interactions between
%  low-complexity data models and Sensing Techniques, Marseille, France,
%  November 21-23, 2018


% <Inputs>
%        X : Input data matrix (m x n)
%        G : Initial G matrix
%        H : Initial H matrix
%        Tmax : Max CPU Time (-0.05) in seconds

% <Default values>
%        MAX_ITER : Maximum number of iterations. Default is 1,000.
%        MIN_ITER : Minimum number of iterations. Default is 10.
%        TOL : Stopping tolerance. Default is 1e-5. If you want to obtain a more accurate solution, decrease TOL and increase MAX_ITER at the same time.

% <Outputs>
%        G : Obtained basis matrix (m x r).
%        H : Obtained coefficients matrix (r x n).
%        T : CPU TIME.
%        RRE: Relative reconstruction error in each iteration

% Note: another file 'stop_rule.m' should be included under the same
% directory as this code.

function [G,H,RRE,T]=RSI_W_NeNMF( X,G,H,r,Tmax)
MinIter=10;
tol=1e-5;
T=zeros(1,301);
RRE=zeros(1,301);

ITER_MAX=500;      % maximum inner iteration number (Default)
ITER_MIN=10;        % minimum inner iteration number (Default)

[L,R]=RSI_compression(X,r);


% Compress left and right
X_L = L * X;
X_R = X * R;


H_comp= H* R;
G_comp = L*G;

HVt=H_comp*X_R';
HHt=H_comp*H_comp';

GtV=G_comp'*X_L;
GtG=G_comp'*G_comp;

GradG=G*HHt-HVt';
GradH=GtG*H-GtV;

init_delta=stop_rule([G',H],[GradG',GradH]);
tolH=max(tol,1e-3)*init_delta;
tolG=tolH; % Stopping tolerance


% Iterative updating
G=G';
k=1;
RRE(k) = nmf_norm_fro( X, G', H);
T(k) =0;
tic
% main loop
while(toc<= Tmax+0.05)
 
    % Optimize H with G fixed
    [H,iterH]=NNLS(H,GtG,GtV,ITER_MIN,ITER_MAX,tolH);
    
    if iterH<=ITER_MIN
        tolH=tolH/10;
    end
    H_comp=H*R;
    HHt=H_comp*H_comp';   HVt=H_comp*X_R';
    % Optimize G with H fixed
    [G,iterG,GradG]=NNLS(G,HHt,HVt,ITER_MIN,ITER_MAX,tolG);
    
    if iterG<=ITER_MIN
        tolG=tolG/10;
    end
    G_comp=G * L';
    GtG=G_comp*G_comp';
    GtV=G_comp*X_L;
    GradH=GtG*H-GtV;
    
    %     HIS.niter=niter+iterH+iterG;
    delta=stop_rule([G,H],[GradG,GradH]);
   
  % Stopping condition
    if (delta<=tol*init_delta && k>=MinIter)
        break;
    end

    
    if toc-(k-1)*0.05>=0.05
        k = k+1;
        RRE(k) = nmf_norm_fro( X, G', H);
        T(k) = toc;
    end
    
end  %end of  loop
G=G';


return;


function [H,iter,Grad]=NNLS(Z,GtG,GtV,iterMin,iterMax,tol)

if ~issparse(GtG)
    L=norm(GtG);	% Lipschitz constant
else
    L=norm(full(GtG));
end
H=Z;    % Initialization
Grad=GtG*Z-GtV;     % Gradient
alpha1=1;

for iter=1:iterMax
    H0=H;
    H=max(Z-Grad/L,0);    % Calculate sequence 'Y'
    alpha2=0.5*(1+sqrt(1+4*alpha1^2));
    Z=H+((alpha1-1)/alpha2)*(H-H0);
    alpha1=alpha2;
    Grad=GtG*Z-GtV;
    
    % Stopping criteria
    if iter>=iterMin
        % Lin's stopping condition
        pgn=stop_rule(Z,Grad);
        if pgn<=tol
            break;
        end
    end
end

Grad=GtG*H-GtV;

return;
function f = nmf_norm_fro(X, G, H)
% Author : F. Yahaya
% Date: 13/04/2018
% Contact: farouk.yahaya@univ-littoral.fr
% Goal: compute a normalized error reconstruction of the mixing matrix V
% "Normalized" means that we divide the squared Frobenius norm of the error
% by the squared Frobenius norm of the matrix V
% Note: To express the error in dB, you have to compute 10*log10(f)
%

f = norm(X - G * H,'fro')^2/norm(X,'fro')^2;

return;

function [ L,R ] = RSI_compression(X,r)
%             Tepper, M., & Sapiro, G. (2016). Compressed nonnegative
%             matrix factorization is fast and accurate. IEEE Transactions
%             on Signal Processing, 64(9), 2269-2283.
%             see: https://arxiv.org/pdf/1505.04650
%             The corresponding code is originally created by the authors
%             Then, it is modified by F. Yahaya.
%             Date: 13/04/2018
%

compressionLevel=20;
[m,n]=size(X);

l = min(n, max(compressionLevel, r + 10));

OmegaL = randn(n,l);

Y = X * OmegaL;

for i=1:4
    
    [Y,~]=qr(Y,0);
    S=X'*Y;
    [Z,~]=qr(S,0);
    Y=X* Z;
end
[L,~]=qr(Y,0);
L=L';

OmegaR = randn(l, m);
Y = OmegaR * X;

for i=1:4
    [Y,~]=qr(Y',0);
    S=X*Y;
    [Z,~]=qr(S,0);
    
    Y=Z'*X;
end
Y=Y';
[R,~] = qr(Y,0);


return
