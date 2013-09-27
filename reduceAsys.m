function Arsys = reduceAsys(Asys)
% Arsys = reduceAsys(Asys)
% ------------------------
%
% Reduces a given Asys cell to a minimal monomial set. Any other monomial
% of Asys can be derived from this mimimal set.
%
% Arsys     =   cell, a polysys cell containing a minimal monomial basis
%               for LT(Id)
%
% CALLS
% -----
%
% Kim Batselier, 2011-06-27
n = size(Asys{1,2},2);

Am = zeros(size(Asys,1),n);
% first put all exponents in a matrix
for i = 1:size(Asys,1)
    Am(i,:) = Asys{i,2}(:)';
end

counter = 1;
while ~isempty(Am)
    
    % first element is always good
    Arsys{counter,1} = 1;
    Arsys{counter,2} = Am(1,:);
    
    verschil = Am-ones(size(Am,1),1)*Arsys{counter,2};
        
    % remove derivatives of baseA{1} from Am
    indices = verschil(:,1) >= 0;
    for i = 2 : n
        indices = indices & (verschil(:,i)>=0);
    end
    Am(indices,:) = [];
    counter = counter + 1;

end