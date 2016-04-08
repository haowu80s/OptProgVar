function [ A, B, R2, objFun, tfocs_out, tfocs_optsOut ] = optimizeAB( X, data_props, C, lam_l1, mu_l2, opts, tol_AB)
%optimizeAB is a driver for  solver_OneSPCC_CVX
%   input:
%       X: cells of data matrix
%       data_props: struct with info of X
%       C: constraints matrix after scaling
%       lam_l1: penalty weight for l1 regularization
%       mu_l2: penalty weight for l2 regularization
%       opts: other options for CVX
%       tol_AB: tolarence for CVX
%   output:
%       A, B, R2: outputs from solver_OneSPCC_CVX. B is the coeff for
%       progress variable before post-processing
%       objFun, tfocs_out, tfocs_optsOut: other outputs for debugging etc.
%%
W = eye(data_props.nPhi);
d = -sqrt(eps)*ones(length(C),1);
[N,~] = size(C);
P = randperm(N,round(N/5));
C_sub = 100*C(P,:);
d_sub = 100*d(P);
%%
[ A, B, R2 ] = solver_OneSPCC_CVX( X, lam_l1, W, C_sub, d_sub, mu_l2, opts, tol_AB );
tfocs_out = [];
tfocs_optsOut = [];
%%
objFun.mse = 0;
for i = 1:length(X)
    objFun.mse = objFun.mse + norm(X{i}-X{i}*B*A{i}', 'fro')^2;
end
objFun.mseNULL = 0;
for i = 1:length(X)
    objFun.mseNULL = objFun.mseNULL + norm(X{i}, 'fro')^2;
end
objFun.l2 = 0.5*mu_l2*norm(B,2)^2;
objFun.l1 = lam_l1*norm(W*B, 1);

objFun.sum = objFun.mse + objFun.l2 + objFun.l1;

end

