function [p d] = elim(polysys,var,varargin)
% [p d] = elim(polysys,var,notsparse)
% -----------------------------------
% Returns a polynomial p which lies in the ideal of the polynomial system
% polysys and from which all variables 'var' have been eliminated.
%
% p         =   row vector, polynomial that lies in the Elimination ideal from
%               which 'var' has been eliminated.
%
% d         =   scalar, degree at which the desired polynomial is found.
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% var       =   vector, contains indices of variables that need to be
%               eliminated.
%
% notsparse =   boolean, if set to 0, then a sparse rank revealing QR
%               is used. If set to 1, a dense SVD is used instead.
%               Default is 0.
% CALLS
% -----
% 
% fite.m, getM.m, getMex.m, updateN.m
%
% Kim Batselier, 2010-12-10, 2011-11 update: uses sparse matrices now
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin == 2
% default behavior is to use sparse data structure and rank revealing QR
    sparse=1; 
else
   sparse=0; 
end

% initialization of outputs
p = [];
d =0;

% number of equations and variables
n_eq = size(polysys,1);
n = size(polysys{1,2},2);
dorig=zeros(1,n_eq);

maxnorm=0;
for i = 1 : n_eq
    polysys{i,1} = polysys{i,1}/norm(polysys{i,1});
    dorig(i)=max(sum(polysys{i,2},2));
    if dorig(i) > d
        d = dorig(i);
    end
    if norm(polysys{i,1}) > maxnorm
        maxnorm = norm(polysys{i,1});
    end
end

mon = getMon(d,n);

% determine the indices where we need nonzero coefficients for the basis
for i = 1 : length(var)    
    ind{i} = find(mon(:,var(i))==0);    
end

temp = ind{1};
for i = 2 : length(var)
    temp = intersect(temp,ind{i});
end

indices = temp;

% initialization Macaulay matrix and orthogonal basis kernel
if sparse
    M=getM(polysys,d,1)';
    % orthogonal basis kernel(M)
    [Q R P]=qr(M,'vector');
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

while isempty(p)
       [Y Sin Z]=svd(full(N(indices,:)'));
       sin=diag(Sin);
       rs=sum(sin > tol);
       
       if (asin(Sin(min(size(Sin)),min(size(Sin)))) < tol) || (rs < size(Sin,2))
			if sparse				
				p=spalloc(1,nchoosek(d+n,n),length(indices));
			else
				p=zeros(1,nchoosek(d+n,n));
			end
			p(1,indices)=Z(:,end)';		
			
       else
			d =d +1;
			clear mon ind temp
			mon = getMon(d,n,d);			
			% determine the indices where we need nonzero coefficients for the basis
			for i = 1 : length(var)        
				ind{i} = find(mon(:,var(i))==0);        
			end
			temp = ind{1};
			for i = 2 : length(var)
				temp = intersect(temp,ind{i});
			end
			indices(length(indices)+1:length(indices)+length(temp)) = 0; 
			for i = 1 : length(temp)
				indices( end-length(temp)+i) = feti(mon(temp(i),:));        
			end    

			if sparse
				N=updateN(N,getMex(polysys,d,d-1,1),1);      
			else
				N=updateN(N,getMex(polysys,d,d-1)); 
			end
       end
end

end
