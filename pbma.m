function [polysys,a,b]=pbma(root,varargin)
% [polysys,a,b]=pbma(root) or pbma(root,mult)
% -------------------------------------------
% Projective Buchberger-Moller algorithm. For a given set of projective roots and their
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
% root      =   matrix, each row corresponds with a projetive root.
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
% getKSB.m, fite.m, getM.m
%
% Kim Batselier, 2014-03

[m,n] = size(root);

if isempty(varargin)
    mult=cell(1,m);
    % no multiplicities
    for i=1:m
        mult{i}=1;
    end
else
    mult=varargin{1};
end

a=[];
b=[2,3];    % 2 pure powers that are connected, hence stop condition not satisfied
polysys=[];
d=0;

while ~canstop(b,n)
    d=d+1;
    % construct kernel K
    % for now, make whole K everytime
    K=[];
    for i=1:m       % for each root
        % determine order of differentiation    
        ddiff=0;
        while nchoosek(ddiff+n,n) < size(mult{i},1)
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
    
    b=[];
    % canonical decomposition
    for i=1:length(indices)
        [U, S,Z]=svd(full(K([b indices(i)],:)'));
        if size(S,2)==1
            S=S(1,1);
        else
            S=diag(S);
        end
        tol=m*S(1)*eps;
        rs=sum(S > tol);
        
        if (S(end) < tol) || (rs < length([b indices(i)]))           
            temp=zeros(1,nchoosek(d+n,n));
            temp([b indices(i)])=Z(:,end);
            temp(abs(temp)<tol)=0;  % remove numerically zero coefficients
            % check whether new polynomial g_i lies in <g_1,..,g_{i-1}>
            if ~isempty(polysys)
               M=getM(polysys,d);               % construct Macaulay matrix
               M=M(:,nchoosek(d-1+n,n)+1:end);  % only need monomials of degree d
               if rank([M;temp(nchoosek(d-1+n,n)+1:end)]) > rank(M)
                   % g_i does not lie in <g_1,..,g_{i-1}>
                   polysys=[polysys;vec2polysys(temp,n)];
                    a=[a indices(i)];
               end
            else
                polysys=[polysys;vec2polysys(temp,n)];
                a=[a indices(i)];
            end
        else
            b=[b indices(i)];
        end
    end    
end

    function satisfied = canstop(b,n)
        satisfied=1;
        
        % determine indices pure powers in B
        exponents=fite(b,n);
        dd=sum(exponents(1,:));
        ppindices=zeros(1,n);
        for ii=1:size(exponents,1) % for each monomial in B            
            if length(find(exponents(ii,:)))==1
                ppindices(find(exponents(ii,:)))=ii;
            end
        end
        ppindices(ppindices==0)=[];  % keep only indices of pure powers in B
        
        for ii=1:length(ppindices)
            % determine connected components for each pure power in B
            connected = b(ppindices(ii));
            index=ppindices(ii);
            
            % going up in indices
            if connected < nchoosek(dd+n,n)
                for jj=[index+1:length(b)]
                    if isconnected(connected,b(jj),n)
                        connected = [connected  b(jj)];
                    end
                end
            end
  
            % going down in indices
            if connected > nchoosek(dd-1+n,n)+1
                for jj=[index-1:-1:1]
                    if isconnected(connected,b(jj),n)
                        connected = [connected  b(jj)];
                    end
                end            
            end
            
            % check divisibility of connected set by x_i
            allexponents=fite(connected,n);
            if ~isempty(find(allexponents(:,ii)==0))
                satisfied=0;
                break;
            end
        end
        
        function isc=isconnected(connected,index,n)
            isc=0;      % is 1 if index is connected
            checkexponents=fite(index,n);
            for iii=1:length(connected)
                bexponents=fite(connected(iii),n);
                % check where there are nonzero exponents in connected(i)
                nzexpindices= find(bexponents >0);
                for jjj=1:length(nzexpindices) % for each nonzero exponent
                    
                    % other exponents
                    otherexponents=1:n;
                    otherexponents(nzexpindices(jjj))=[];
                                        
                    % construct multiply and divide exponents
                    mdexponents=zeros(n-1,n);
                    mdexponents(:,nzexpindices(jjj))=-1;
                    for k=1:n-1
                        mdexponents(k,otherexponents(k))=1;
                    end
                    
                    % check reach of bexponents
                    for k=1:n-1
                        if sum(bexponents+mdexponents(k,:) == checkexponents)==n
                           % the monomial is connected to pure power
                           isc=1;
                           break;
                        end                        
                    end                    
                end
            end
            
        end
        
    end

end