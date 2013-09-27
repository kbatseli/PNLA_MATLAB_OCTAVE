function sol = makeRoot(d,root)
% sol = makeRoot(d,root)
% ----------------------
%
% Evaluates a multivariate polynomial base vector.
%
% sol   =   column vector, contains the evaluated polynomial base vector
%
% d     =   scalar, degree of the multivariate polynomial base vector
%
% root  =   row vector, used to evaluate differentiated base with
%
% CALLS
% -----
%
% getMon.m
%
% Kim Batselier, 2010-01

n = size(root,2);

monBase = getMon(d,n);

l = length(monBase);

for j = 1 : size(root,1)
    temp = zeros(l,n);
    for i = 1 : n
        temp(:,i) = (root(j,i)*ones(l,1)).^monBase(:,i);
    end
    
    sol(:,j) = prod(temp,2);
end

end