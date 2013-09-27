function D = getD(polysys,p,varargin)
% D = getD(polysys,p,sparseM)
% ----------------------------
% 
% This function makes the divisor matrix of given set of polynomial equations
% polysys and polynomial p of degree d.
%
% D         =   Divisor matrix D,
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% p         =   vector, contains coefficients of polynomial p
%
% sparseM   =   boolean, optional: set to 1 if Qpd needs to be sparse
%               represented, default: 0
%
% CALLS
% -----
%
% getMon.m, feti.m, getMonBase.m
%
% Kim Batselier, 2011-03-15,16
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

if nargin < 3
    sparseM = 0;
else
    sparseM = varargin{1};
end

% ensure p is a row vector
p = p(:)';

% number of variables
n = size(polysys{1,2},2);
% number of equations
n_eq = size(polysys,1);
% this vector will contain the degree of each equation
dorig = zeros(n_eq,1);
monorig = zeros(n_eq,n);
shiftind = zeros(n_eq,1);

% need to know leading term of p
maxcol = find(p,1,'last');

% degree of p
d = sum(fite(maxcol,n));

do =0;

% initialize number of rows of M
r= n_eq;
ncoef=zeros(1,n_eq);

for i = 1 : n_eq
   
    dorig(i,1) = max(sum(polysys{i,2},2));
    if max(sum(polysys{i,2},2)) > do
        do = max(sum(polysys{i,2},2));
    end
    for j = 1 : size(polysys{i,2},1)
        if feti(polysys{i,2}(j,:)) > feti(monorig(i,:))
            monorig(i,:) = polysys{i,2}(j,:);
        end
    end
    
    % determine whether we can shift or not
    if feti(monorig(i,:)+[1 zeros(1,n-1)]) > maxcol
        % shfiting with merely x_1 already puts us over the maximum column
        shiftind(i) = 0;
    else
        % now determine how many shifts this polynomial needs
        addBase = getMon(d-dorig(i),n);
        addBase(1,:) = [];
        for j = 1 : size(addBase,1)
            if feti(monorig(i,:)+addBase(j,:)) <= maxcol
                shiftind(i) = shiftind(i)+1;
            else
                break
            end
        end        
    end
   
    r = r + shiftind(i);
     ncoef(i) = length(polysys{i,1}) *(1 + shiftind(i));
end

% check degree p with max degree of polysys
if d < do
    error('Degree of given polynomial p is smaller than polysys.')
end

% number of columns of D
c = maxcol;

% first allocate memory for the divisor matrix to speed up things
if sparseM
    D = sparse([],[],[],r,c,sum(ncoef));
else
    D = zeros(r,c);
end

rowcounter = 1;

shiftd = d-dorig;   % max degree of shift monomial


for i =1: n_eq    % for each equation
    
    % determine the monomials we additionally need to multiply with
    addBase = getMon(shiftd(i,1),n);
    
    for j = 1 : size(addBase(1:shiftind(i)+1,:),1)     % for each shift
        col=zeros(1,size(polysys{i,2},1));
        
        for k = 1 : size(polysys{i,2},1) % for each monomial in the equation
            
            col(k) = feti(addBase(j,:)+polysys{i,2}(k,:));            
        end        
        
        D(rowcounter,col) = polysys{i,1};
        rowcounter = rowcounter + 1;
    end
end


end
