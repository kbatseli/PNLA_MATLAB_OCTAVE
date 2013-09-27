function [di n] = getDorig(polysys)
% [di n] = getDorig(polysys)
% -----------------------------------
%
% Returns the degrees of each polynomial in polysys.
%
% di        =   vector, each entry is the degree of the corresponding
%               polynomial in polysys,
%
% n         =   scalar, number of variables,
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               given polynomial system.
%
% CALLS
% -----
%
% Kim Batselier, 2011-05
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


n = size(polysys{1,2},2);
n_eq = size(polysys,1);
di = zeros(n_eq,1);

for i = 1 : n_eq    
    di(i,1) = max(sum(polysys{i,2},2));
end

end
