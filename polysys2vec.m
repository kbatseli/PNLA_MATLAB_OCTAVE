function y = polysys2vec(polysys,d,varargin)
% y = polysys2vec(polysys,d,sparsep)
% ----------------------------------
% 
% Takes every equation in the polysys cell and converts it to a vector in a
% certain monomial basis.
%
% y         =   matrix, every row vector corresponds with an equation from
%               the polysys cell.
%
% polysys   =   cell, same structure as the polysys cells for getM
%
% d         =   scalar, maximal degree of the monomial basis which is
%               used. If d > max(d_i) then zeros are appended.
%
% sparsep   =   boolean, set to 1 if result needs to be sparse, default = 0
%
% CALLS
% -----
%
% getMon.m, feti.m
%
% Kim Batselier, 2009-11-13, update 2012: added 'degree' argument

n = size(polysys{1,2},2);
n_eq = size(polysys,1);
di = zeros(n_eq,1);
sparsep = 0;
maxl = 0;   

if nargin == 3
    sparsep = varargin{1};     
end

for i = 1 : n_eq    
    di(i,1) = max(sum(polysys{i,2},2));
    if length(polysys{i,1}) > maxl
        maxl = length(polysys{i,1});
    end
end

if d < max(di)
    error(['Polynomial(s) ' num2str(find(di <= d)) ' have a degree which is larger than d']);
end

if sparsep
    y = sparse([],[],[],size(polysys,1),nchoosek(d+n,n),maxl);
else
    y = zeros(n_eq,nchoosek(d+n,n));
end

rowcounter = 1;
% first we add shifted versions of all equations with degree < do so that
% they all have degree do
for i =1: n_eq    % for each equation
    col = zeros(1,size(polysys{i,2},1));
    for k = 1 : size(polysys{i,2},1) % for each monomial in the equation
        
        col(k) = feti(polysys{i,2}(k,:));            
    end
    
    y(rowcounter,col) = polysys{i,1};
    rowcounter = rowcounter + 1;

end

end
