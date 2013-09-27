function exponent = fite(index,n)
% exponent = fite(index,n)
% ------------------------
% 
% Converts a vector of indices with respect to the degree negative lex
% monomial ordering into a matrix of exponents.
%
% index     =   vector, each entry corresponds with an index of a
%               particular monomial, 
%
% n         =   scalar, number of variables of each of the monomials,
%
% exponent  =   matrix, each row exponent(i,:) is the exponent of the
%               monomial with index index(i).
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

exponent=zeros(length(index),n);

for j=1:length(index)
    
    % first asses the degree
    d=0;
    while nchoosek(d+n,n) < index(j)
        d=d+1;
    end
    
%     if d==0
%         break
%     else
     if d >0   
        index(j) = index(j)-nchoosek(d-1+n,n);
        
        for i=1:n-1
            k=0;
            while nchoosek(k+n-i-1,n-i-1) < index(j)
                index(j) = index(j) - nchoosek(k+n-i-1,n-i-1);
                k=k+1;
            end
            
            exponent(j,i)=d-sum(exponent(j,1:i-1))-k;
        end
        
        exponent(j,end) = d-sum(exponent(j,1:end-1));
        
    end
end

end
