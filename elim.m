function [p d] = elim(polysys,var)
% [p d] = elim(polysys,var)
% -------------------------
% Returns a polynomial p which lies in the ideal of the polynomial system
% polysys and from which all variables 'var' have been eliminated.
%
% p         =   row vector, polynomial that lies in the Elimination ideal from
%               which 'var' has been eliminated
%
% d         =   scalar, degree of the polynomial from which variables were
%               eliminated
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% var       =   vector, contains indices of variables that need to be
%               eliminated
%
% CALLS
% -----
% 
% getM.m, getMon.m, getMex.m
%
% Kim Batselier, 2010-12-10, 2011-11 update: uses sparse matrices now
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

% initialization of outputs
p = [];
d =0;

% number of equations and variables
n_eq = size(polysys,1);
n = size(polysys{1,2},2);
dorig=zeros(1,n_eq);

maxnorm=0;
for i = 1 : n_eq
    polysys{i,1} = polysys{i,1}/norm(polysys{i,1});
    dorig(i)=max(sum(polysys{i,2},2));
    if dorig(i) > d
        d = dorig(i);
    end
    if norm(polysys{i,1}) > maxnorm
        maxnorm = norm(polysys{i,1});
    end
end

mon = getMon(d-1,n);

% determine the indices of the nonzero coefficients of the monomial
% basis of E(d)
for i = 1 : length(var)    
    ind{i} = find(mon(:,var(i))==0);    
end
temp = ind{1};
for i = 2 : length(var)
    temp = intersect(temp,ind{i});
end
indices = temp;


[M nztemp] = getM(polysys,d-1,1);

while isempty(p)
    
    clear mon ind temp
    mon = getMon(d,n,d);
    
    % determine the indices of the nonzero coefficients of the monomial
    % basis of E(d)
    for i = 1 : length(var)        
        ind{i} = find(mon(:,var(i))==0);        
    end    
    temp = ind{1};
    for i = 2 : length(var)
        temp = intersect(temp,ind{i});
    end
    
    indices(length(indices)+1:length(indices)+length(temp)) = 0; 
    for i = 1 : length(temp)
        indices( end-length(temp)+i) = fetr(mon(temp(i),:));        
    end    
    
    c = nchoosek(d+n,n);
    r= n_eq;
    for i = 1 : n_eq
        r = r + nchoosek(d-dorig(i)+n,n)-1;
    end
    
    tau = 20*(r+c)*maxnorm*eps;
    
    % update M
    temp = M;
    [Mex nzmex] = getMex(polysys,d,d-1,1);
    M = spalloc( r,c ,nztemp+nzmex);
    M(1:size(temp,1),1:size(temp,2)) = temp;
    M(size(temp,1)+1:end,:) = Mex;
    clear Mex

    [Q, R,P] = qr(M',0);
    rankM = length(find(diag(R)));    
    U = Q(:,1:rankM); % basis for row space
        
    % sine-based angle calculation
    test=-U*U(indices,:)';
    for i =1:length(indices)
        test(indices(i),i) = 1+test(indices(i),i);
    end
    [~, C Z] = svd(full(test));
    s=diag(C);

    if asin(s(end))<tau
        % sine-based principal vector
        p=spalloc(1,size(U,1),length(indices));
        p(indices)=Z(:,end);
        break        
    else
        d = d+1;
    end
    
end

end
