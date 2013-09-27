function [ang p] = idealMember(psys,polysys,d,varargin)
% [ang p] = idealMember(psys,polysys,d,tol)
% -----------------------------------------
% Ideal membership checking via principal angles. Returns principal angle 
% between M(d) and monomial basis for given polynomial psys. Computes sines
% of principal angles instead of cosines.
%
% ang	 	= 	scalar, principal angle between row space M(d) and given
%               polynomial
%
% p	    	=   vector, coefficient vector of polynomial that lies in row(M(d))
%
% psys  	=   cell, polysys cell for polynomial you want to check
%
% polysys 	=	cell, polysys cell for generators of ideal
%
% d		    =   scalar, degree of row(M(d))
%
% tol       =   scalar, tolerance: optional. Default =
%               max(size(M))*eps(sigma1) with sigma1 = largest singular value
%
% CALLS
% -----
% 
% getM.m, feti.m
%
% Kim Batselier, 2012-05-28

M = getM(polysys,d);
[A S V] = svd(M,'econ');
sig=diag(S);
if nargin == 3
    tol = max(size(M))*eps(sig(1));
else
    tol = varargin{1};
end

% need monomial basis for psys
indices = zeros(1,length(psys{1,1}));
for i = 1:length(psys{1,1})
    indices(i) = feti(psys{1,2}(i,:));
end
%E = speye(fetr(LT),size(M,2));
indices = sort(indices);
p=zeros(1,size(M,2));

r = sum(sig>tol);
U = V(:,1:r);
clear A S V

test=-U*U(indices,:)';
for j=1:length(indices)
    test(indices(j),j) = 1+test(indices(j),j);
end
% compute the sines
[Y S Z] = svd(test);
s=diag(S);

ang = asin(s(end));

if  ang < tol		
 p(1,indices)=Z(:,end)';
end

end
