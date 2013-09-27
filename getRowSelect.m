function [Sa Sb] = getRowSelect(d,n,a,b)
% [Sa Sb] = getRowSelect(d,n,a,b)
% -------------------------------
% Returns row selection matrices for shifts a,b of homogeneous n-variate Vandermonde vectors of degree d.
%
% Sa	        =   vector, contains indices of left-hand side row selection, multiplication with x_a
%
% Sb	        =   vector, contains indices of right-hand side row selection, multiplication with x_b
%
% d		=   scalar, degree at which homogeneous Vandermonde vector is made
%
% n             =   scalar, number of variables
%
% a 		=   scalar, index of multiplication monomial on left-hand side, 0 <= a <=n
%
% b 		=   scalar, index of multiplication monomial on right-hand side, 0 <= b <=n
%
% CALLS
% -----
% 
% getMonBase.m
%
% Kim Batselier, 2013-09


if a<0 | a > n
    error(['Error: shift variable should be scalar between ' num2str(0) ' and ' num2str(n)])
end

if b<0 | b > n
    error(['Error: shift variable should be scalar between ' num2str(0) ' and ' num2str(n)])
end

% monomial basis
mons=getMonBase(d,n+1);

% left-hand side row selection matrix depends on right-hand side shift
% monomial b
sa=find(mons(:,b+1)>0);
Sa = zeros(length(sa),nchoosek(d+n,n));
for i=1:length(sa)
    Sa(i,sa(i))=1;
end

% right-hand side row selection matrix depends on left-hand side shift
% monomial b
Sa=find(mons(:,b+1)>0)';
Sb=zeros(1,length(Sa));
for i=1:length(Sa)
    indexfound=0;
    index=1;
    while ~indexfound
        if sum(mons(index,:)==mons(sa(i),:)+[zeros(1,a) 1 zeros(1,n-a)]-[zeros(1,b) 1 zeros(1,n-b)])==n+1
            Sb(i)=index;
            indexfound=1;
        else
            index=index+1;
        end
    end            
end

end
