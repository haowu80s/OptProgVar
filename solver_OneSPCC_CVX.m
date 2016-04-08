function [ A, B, R2 ] = solver_OneSPCC_CVX( X, lambda, W, C, d, mu, opts, tol_AB )
% solver_OneSPCC_CVX Minimize residual subject to l1-norm constraints,
%  and l2-norm penalty.
% [ A, B, out, opts ] = solver_SPCC( X, k, lambda, mu, B0, opts, varargin )
%    Solves the SPCC modified from the SPCA formulation the Zou & Tibshirani formulation,
%        minimize sum_theta sum_i norm(x_i-A(theta)B^Tx_i,2)^2+ sum_j lambda_j * norm( W*beta_j, 1 ) +
%        0.5*mu*norm((beta_j-beta_j0),2)^2
%
% See also solver_L1RLS.m, solver_sBPDN.m, solver_EN from CVX

% define residual norm
spca_out.lp='fro';
% check size of X
m = length(X);
[~,p] = size(X{1});
n = zeros(m,1);
for i=1:m
    [n(i), ~] = size(X{i});
end
% check input parameter
if nargin < 6 || isempty(mu), mu = sqrt(eps); end
if nargin < 7, opts = []; end
%% define solution format
A{m} = [];
S{m} = [];
B = zeros(p,1);
%%
% define residual norm
if (~isfield( opts, 'spca_out.lp' ))
    spca_out.lp=Inf;
end
% set initial values
for i=1:m
    [~,~,V] = svd(X{i},'econ');
    A{i} = V(:,1);
    S{i} = X{i}'*X{i};
end
clear V;
% set X2
X2=zeros(sum(n),p);
n_cum = cumsum(n);
n_cum = [0;n_cum];
for i=1:m
    X2(n_cum(i)+1:n_cum(i+1),:) = X{i};
end
%%
% initialze output

disp('     it      step-A      step-B      step-F');
disp('     --------------------------------------');
% iterate betwen A and B
spca_out.nIt = 0;
spca_out.res_B = 1;
spca_out.res_A = 1;
objFun.mse = 0;
%%
while (true)
    spca_out.nIt = spca_out.nIt+1;
    B_old = B;
    
    Y2=zeros(sum(n),1);
    for i=1:m;
        Y2(n_cum(i)+1:n_cum(i+1)) = X{i}*A{i};
    end
    %%
    B = solver_LICEN_CVX( X2, Y2, C, d, mu, lambda, W );
    if norm(B,2)<0.9
        error('Penalty is too large')
    end
    %% find A given B
    A_old = A;
    for i=1:m
        %[U,~,V] = svd(S{i}*B, 'econ');
        %A(:,:,i) = U*V';
        A_tmp = S{i}*B;
        if norm(A_tmp,2)<norm(B,2)*sqrt(eps)
            A{i} = ones(p,1);
            warning('norm(A_tmp,2)<norm(B,2)*sqrt(eps)');
        end
        A{i} = A_tmp/(norm(A_tmp,2)+sqrt(eps));
    end
    % res_AB
    spca_out.res_B = norm(B-B_old,spca_out.lp)./ norm(B,spca_out.lp);
    res_A_tmp = 0;
    for i=1:m
        res_A_tmp = res_A_tmp + sqrt(sum(sum((A{i}-A_old{i}).^2)))./ sqrt(sum(sum((A{i}).^2)));
    end
    spca_out.res_A = res_A_tmp;
    % res_F
    objFun.mse_old = objFun.mse;
    objFun.mse = 0;
    for i = 1:length(X)
        objFun.mse = objFun.mse + norm(X{i}-X{i}*B*A{i}', 'fro')^2;
    end
    objFun.l2 = 0.5*mu*norm(B,2)^2;
    objFun.l1 = lambda*norm(W*B, 1);
    objFun.sum = objFun.mse + objFun.l2 + objFun.l1;
    res_F = abs(objFun.mse-objFun.mse_old)/abs(objFun.mse);
    disp([num2str(spca_out.nIt), ' ', num2str(spca_out.res_A,'%1.2e'), '  ',...
        num2str(spca_out.res_B,'%1.2e'), '  ', num2str(res_F,'%1.2e')]);
    if (spca_out.res_A<=tol_AB && spca_out.res_B<=tol_AB ) || (res_F < tol_AB)
        break
    end
    
end
% Adjusted Totoal Variance
R2 = diag(qr(X2*B,0)).^2;
% Total Varriance
spca_out.TV=trace(X2'*X2);

