function [puni d sintheta tol] = punivar(polysys,var,varargin)
% [puni d sintheta tol] = punivar(polysys,var,notsparse)
% ------------------------------------------------------
% Returns the univariate polynomial p(x_i) obtained by the elimination of
% all other variables from the polynomial system polysys.
%
% puni          =   row vector, coefficient vector of the univariate
%                   polynomial, degree increases from left to right
%
% d 			=	scalar, degree at which the desired polynomial is found.
%
% sintheta  	=	scalar, sine of the smallest principal angle between the
%					row space of M(d) and the eliminiation space.
%
% tol 			=	scalar, numerical tolerance used for the rank-test.
%
% polysys       =   cell containing coefficients and monomials exponents of the
%                   set of polynomial equations.
%
% var           =   scalar, the index i of the variable x_i of the
%                   required univariate polynomial
%
% notsparse     =   boolean, if set to 0, then a sparse rank revealing QR
%                   is used. If set to 1, a dense SVD is used instead.
%                   Default is 0.
%
% CALLS
% -----
% 
% getD0.m, getM.m, updateN.m
%
% Kim Batselier, 2013

if nargin == 2
% default behavior is to use sparse data structure and rank revealing QR
    sparse=1; 
else
   sparse=0; 
end

% number of variables
n=size(polysys{1,2},2);
d=getD0(polysys);
sintheta=0;

for i=1:size(polysys,1)
    polysys{i,1}=polysys{i,1}/norm(polysys{i,1});
end

% initialization Macaulay matrix and orthogonal basis kernel
if sparse
    M=getM(polysys,d,1)';
    % orthogonal basis kernel(M)
    [Q R P]=qr(M);
    r=nnz(diag(R));
    N=Q(:,r+1:end);
    tol=20*sum(size(M))*eps;
else
    M=getM(polysys,d)';
    [U S V]=svd(M);
    s=diag(S);
    tol=max(size(M))*eps(s(1));
    r=sum(s > tol );
    N=U(:,r+1:end);
end

% initialization monomial basis for univariate polynomial
indices=ones(1,d+1);
for i=1:d
    indices(i+1)=feti([zeros(1,var-1) i zeros(1,n-var)]);
end

puni=[];

while isempty(puni)
       [Y Sin Z]=svd(full(N(indices,:)'));
       
       if (asin(Sin(min(size(Sin)),min(size(Sin)))) < tol) || (size(Sin,1) < size(Sin,2))
			puni=Z(:,end)';
			if size(Sin,1) >= size(Sin,2)
				sintheta=Sin(min(size(Sin)),min(size(Sin)));
			end
       else
			d =d +1;
			indices=[indices feti([zeros(1,var-1) d zeros(1,n-var)])];
			if sparse
				N=updateN(N,getMex(polysys,d,d-1,1),1);      
			else
				N=updateN(N,getMex(polysys,d,d-1)); 
			end
       end
end


end
