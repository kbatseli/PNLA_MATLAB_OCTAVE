function polysys = vec2polysys(vec,n)
% polysys = vec2polysys(vec,n)
% ----------------------------
% 
% Converts a vector with coefficients of a multivariate polynomial into a
% polysys cell.
%
% polysys   =   cell, same structure as the polysys cells for getM
%
% vec       =   matrix, each row is a coefficient vector of a multivariate
%               polynomial
%
% CALLS
% -----
%
% fite.m
%
% Kim Batselier, 2009-11-13, 2011-08-03, Rewritten: removed ordering
% argument, vec input can now be a matrix, each row of the matrix is then
% interpreted as the coefficient vector of a polynomial. Input needs to be
% a row vector!
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

if vec==0
    polysys = cell(1,2);
    polysys{1,1} = 0; polysys{1,2} = zeros(1,n);
else
    
polysys = cell(size(vec,1),2);

for i = 1 : size(vec,1) % each row of vec is a coefficient vector
    coefI = find(vec(i,:));
    polysys{i,1} = vec(i,coefI);
    for j = 1 : length(coefI)
        polysys{i,2} = [polysys{i,2} ; fite(coefI(j),n) ];
    end
end


end

end
