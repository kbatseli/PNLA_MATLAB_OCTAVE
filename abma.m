function [polysys,a,b]=abma(root,varargin)
% [polysys,a,b]=abma(root) or abma(root,mult)
% -------------------------------------------
% Affine Buchberger-Moller algorithm. For a given set of roots and their
% multiplicity structures, this algorithm returns a mimimal reduced Grobner
% basis.
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               reduced Grobner basis.
%
% a         =   vector, indices of monomials which are leading monomials
%               that can be reached in C_d^n.
%
% b         =   vector, indices of standard monomials.
%
% root      =   matrix, each row corresponds with an affine root.
%
% mult      =   cell, mult{i} contains a matrix T that defines the
%               multiplicit structure of root(i,:). The matrix T is such
%               that the canonical kernel K = D*T, with D the matrix of
%               differential functionals obtained via getKSB.m. Optional,
%               default is no multiplicities.
%
% CALLS
% -----
% 
% getKSB.m, fite.m
%
% Kim Batselier, 2014-03

a=[];
b=1;
polysys=[];

[m,n] = size(root);

if isempty(varargin)
    % no multiplicities
    for i=1:m
        mult{i}=1;
    end
end

stop=0;
d=0;

while ~stop
    d=d+1;
    % construct kernel K
    % for now, make whole K everytime
    K=[];
    for i=1:m       % for each root
        % determine order of differentiation    
        ddiff=0;
        while nchoosek(d+n,n) < size(mult{i},1)
            ddiff=ddiff+1;
        end
        
        D=getKSB(d,ddiff,root(i,:));
        
        K=[K D*mult{i}];
    end
        
    indices = nchoosek(d-1+n,n)+1:nchoosek(d+n,n); % indices of all monomials of degree d
    
    % remove multiples of A from indices
    for i=1:length(indices) % for each new monomial
        for j=1:length(a)   % check whether it is multiple of a(j)
            if  sum((fite(indices(i),n)-fite(a(j),n)) >= 0) == n
                % we found a multiple
                indices(i) = 0;
            end
        end        
    end
    indices(indices==0)=[];
    if isempty(indices)
        stop=1;
    else
        % canonical decomposition
        for i=1:length(indices)
            [~, S Z]=svd(full(K([b indices(i)],:)'));
            if size(S,2)==1
                S=S(1,1);
            else
                S=diag(S);
            end
            tol=m*S(1)*eps;
            rs=sum(S > tol);
            
            if (S(end) < tol) || (rs < length([b indices(i)]))
                a=[a indices(i)];
                temp=zeros(1,nchoosek(d+n,n));
                temp([b indices(i)])=Z(:,end);
                temp(abs(temp)<tol)=0;  % remove numerically zero coefficients
                polysys=[polysys;vec2polysys(temp,n)];
            else
                b=[b indices(i)];
            end
        end
    end
end

end