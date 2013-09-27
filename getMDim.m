function [r c] = getMDim(polysys,d)
% [r c] = getMDim(polysys,d)
% --------------------------
% Returns the dimensions of the Macaulay matrix M of a polynomial system
% for a given degree d.
%
% r         =   scalar, number of rows of M
%
% c         =   scalar, number of columns of M
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% d         =   scalar, maximum degree of M
%
% CALLS
% -----
% 
% Kim Batselier, 2009-10-13
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


% number of unknowns
n = size(polysys{1,2},2);

% number of equations
n_eq = size(polysys,1);

% get the full base of monomials for the original set of equations, need
% maximum degree present in the set of equations
dmin =0;
r= n_eq;                                       % number of rows of M
for i = 1 : n_eq
    dorig(i,1) = max(sum(polysys{i,2},2));   
    if max(sum(polysys{i,2},2)) > dmin
        dmin = max(sum(polysys{i,2},2));
    end
end

% check given degree argument
if d < dmin
%     warning('You have provided a degree smaller than the maximum total degree of the given set of equations. Setting the degree to the maximum total degree in the set of equations...')
    % need to consider only the polynomials of degree d or smaller
    indices = find(dorig <= d);
    if ~isempty(indices)
        for i = 1 : length(indices)
            temp{i,1} = polysys{indices(i),1};
            temp{i,2} = polysys{indices(i),2};
        end
        clear polysys
        polysys = temp;
        dorig = dorig(indices,1);
        n_eq = length(indices);
    else
        r=0;
        c=0;
        return
    end
end

% initialize number of rows of M
r= n_eq;
for i = 1 : n_eq
    r = r + nchoosek(d-dorig(i)+n,n)-1;
end

% number of columns of M
c = nchoosek(d+n,n);

end

