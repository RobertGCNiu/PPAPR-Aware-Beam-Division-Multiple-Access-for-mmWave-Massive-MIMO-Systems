function x  = cvxToOptimize( He, fbb, lamda )
n = size(fbb,1);
cvx_begin
        variable x(n)
        minimize(norm(He*x-He*fbb,'fro'))
        subject to
                abs(x)<=lamda;
cvx_end
end


%A = randn(5,5);
% b=randn(5,5)+1;
% cvx_begin
%         variable A(5,5)
%         minimize(norm(A-b));
%         subject to
%                 abs(A)<=1;
%  cvx_end
%  A