function index=feti(exponent)
% index=feti(exponent)
% --------------------
% 
% Converts a matrix of exponents of n-variate monomials to their
% corresponding indices with respect to the degree negative lexicographic
% monomial ordering.
%
% index     =   vector, each entry corresponds 
%
% exponent  =   vector, exponent of a monomial
%
% CALLS
% -----
%
% Kim Batselier 2013-07
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

% k is the number of monomials, n the number of variables
[k,n]=size(exponent);
d=sum(exponent,2);

index=zeros(1,k);

for i=1:k
    if sum(exponent(i,:))==0
        index(1,i)=1;
    else
            index(1,i)=nchoosek(d(i)-1+n,n);
            for j=1:n-2
                if d(i)-sum(exponent(i,1:j))-1+n-j >= n-j
                    index(1,i) = index(1,i) + nchoosek(d(i)-sum(exponent(i,1:j))-1+n-j,n-j);
                end
            end
            index(1,i)=index(1,i) + exponent(i,end)+1;
    end
end
end