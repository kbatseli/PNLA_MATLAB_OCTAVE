function [a b R varargout] = rcandecomp(polysys,d,varargin)
% [a b R] = rcandecomp(polysys,d,varargin)
% ----------------------------------------
% Computes the reduced canonical decomposition for a given polynomial
% system polysys at degree d. Also returns for each leading monomial in the
% row space the corresponding polynomial in the ideal. If the degree d is
% large enough, this will correspond with a reduced Groebner Basis.
%
% a         =   vector, indices of reduced monomials which are leading monomials
%               that can be reached in C_d^n.
%
% b         =   vector, indices of monomials which lie in the complement of
%               vector space spanned by monomials of a = affine normal set.
%
% R         =   sparse matrix, each row i corresponds with a polynomial which has as
%               leading monomial the corresponding monomial of a(i).
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% d         =   scalar, degree for which the canonical decomposition is
%               computed.
%
% tol       =   scalar, optional tolerance for checking numerical zeros.
%               Default: sum(size(N))*eps.
%
% CALLS
% -----
% 
% getM.m, updateN.m, getMex.m,
%
% Kim Batselier, 2013-07
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

% initialize
n=size(polysys{1,2},2);
d0=getD0(polysys);
N=null(getM(polysys,d0));
angles=zeros(1,size(N,1));

% M=getM(polysys,d0,1);
% [Q R P]=qr(M','vector');
% r=nnz(diag(R));
% c(d)=size(M,2)-r;
% N=Q(:,r+1:end); 
% clear Q R P

% recursively update orthogonal basis for null space
for i=d0+1:d
    N=updateN(N,getMex(polysys,i,i-1,1),1);
end
c=size(N,2); % corank

if isempty(varargin)
    tol=sum(size(N))*eps;
else
    tol=varargin{1};
end

% first check whether 1 is in the ideal
[Y,Sin,Z]=svd(full(N(1,:)'));
if size(Sin,2)==1
    sin=Sin(1,1);
else
    sin=diag(Sin);
end
rs=sum(sin > tol);
  
if (asin(Sin(min(size(Sin)),min(size(Sin)))) < tol) || (rs < size(Sin,2))
    % 1 is in the ideal, we can quit here
    a=1;
    b=[];
    R=[];
else
    checki = 1:size(N,1); %contains indices of monomials that need to be checked for linear independence
    a=zeros(1,nchoosek(d+n,n)-1);
    acounter=1;

    R=spalloc(nchoosek(d+n,n),nchoosek(d+n,n),nchoosek(d+n,n)*(c+1));
    rcounter=1;
    
    b=zeros(1,nchoosek(d+n,n));
    b(1)=1;
    bcounter=2;
    %         b=[b checki(counter)];
end

% counter=1;
% while counter <= length(checki)
for counter=2:length(checki)
    if counter > length(checki)
        break
    end
    [Y,Sin,Z]=svd(full(N([b(1:bcounter-1) checki(counter)],:)'));
    if size(Sin,2)==1
        sin=Sin(1,1);
    else
        sin=diag(Sin);
    end
    rs=sum(sin > tol);
    
    
    if (asin(Sin(min(size(Sin)),min(size(Sin)))) < tol) || (rs < size(Sin,2))
%         a=[a checki(counter);]
        a(acounter) = checki(counter);
        angles(acounter)=asin(Sin(min(size(Sin)),min(size(Sin))));
        acounter=acounter+1;
        % remove all monomial multiples from checki(counter)
        di=sum(fite(checki(counter),n));
        if di<d
            multiplei=2:nchoosek(d-di+n,n); % indices of monomial multiples
            for j=1:length(multiplei)
                checki(checki==feti(fite(checki(counter),n)+fite(multiplei(j),n)))=[];
%                 length(checki)
            end
        end
        R(rcounter,[b(1:bcounter-1) checki(counter)])=Z(:,end);
        rcounter=rcounter+1;
    else
        b(bcounter)=checki(counter);
        bcounter=bcounter+1;
%         b=[b checki(counter)];
    end
%     counter = counter + 1;
end
b=b(1:bcounter-1);
a=a(1:acounter-1);
angles=angles(1:acounter-1);
R=R(1:acounter-1,:);
varargout{1}=angles;
end


