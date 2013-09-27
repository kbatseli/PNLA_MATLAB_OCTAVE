function d = getD0(polysys)
% d = getD0(polysys)
% ------------------
% Returns the maximum total degree of the given polynomial system
%
% d         =   scalar, maximum total degree of the given polynomial
%               system,
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               given polynomial system.
%
%
% CALLS
% -----
%
% Kim Batselier, 2010-04-21
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


% number of variables
n = size(polysys{1,2},2);

% number of equations
n_eq = size(polysys,1);

% get the full base of monomials for the original set of equations, need
% maximum degree present in the set of equations
d =0;

for i = 1 : n_eq
    if max(sum(polysys{i,2},2)) > d
        d = max(sum(polysys{i,2},2));
    end
end

end