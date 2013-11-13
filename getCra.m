function rows = getCra(Asys,dmax)
% rows = getCra(Asys)
% -------------------
% 
% Returns the indices of the rows of the numerical kernel to construct the
% eigenvalue problem of only the affine roots.
%
% Asys      =   cell, containing coefficients and monomials exponents of the
%               set of monomial system.
% CALLS
% -----
%
% getMon.m, feti.m
%
% Kim Batselier, 2011-10

n_eq = size(Asys,1);
n = size(Asys{1,2},2);
% dorig = zeros(1,n_eq);
% d = 0;
% for i = 1 : n_eq
%     dorig(1,i) = sum(Asys{i,2},2);
% %     if max(sum(Asys{i,2},2)) > d
% %         d = max(sum(Asys{i,2},2));
% %     end
% end
% d = sum(dorig);

monBase = getMon(dmax,n);

for i = 1 : n_eq
    % remove all derivatives of Asys from monBase    
    verschil = monBase-ones(size(monBase,1),1)*Asys{i,2};
        
    indices = verschil(:,1) >= 0;
    for i = 2 : n
        indices = indices & (verschil(:,i)>=0);
    end
    monBase(indices,:) = [];
end

rows = zeros(1,size(monBase,1));
for i = 1 : length(rows)
    rows(i) = feti(monBase(i,:));
end


end
