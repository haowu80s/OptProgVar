function [ x ] = solver_LICEN_CVX( A, b, C, d, mu, lambda, W )
%%
[~,n] = size(A);
%%
cvx_begin quiet
cvx_precision default
   variable x(n)
   dual variable y
   minimize( 0.5*(A*x-b)'*(A*x-b)+ 0.5*mu*x'*x +  lambda*norm(W*x,1)); %/scale + 0.00001*x'*x + 0.0*norm(W*x,1)
    subject to
         y :  C*x >= d;
         x(1) == 0;
cvx_end

%%
% subplot(3,1,1)
% plot(C_sub*x,'.')
% subplot(3,1,2)
% plot(C*x,'.')
% x(abs(x)<1e-2) = 0;
% subplot(3,1,3)
% bar(abs(x))
% sum(x~=0)
%%
end

