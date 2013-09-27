function [hasAll components indices] = hasAllPureComponents(Asys,varargin)
% [hasAll components indices] = hasAllPureComponents(Asys) or hasAllPureComponents(a,n)
% -------------------------------------------------------------------------------------
% Checks whether Asys contains the pure components of every indeterminate for a certain exponent.
%
% hasAll        =   boolean, returns 1 if pure component for each indeterminate
%                   is present in Asys
%
% components    =   row vector, contains for each x_i the exponent of its
%                   pure component in Asys
%
% indices       =   row vector, contains for each x_i the index in Asys
%                   where its corresponding pure component is located.
%
% Asys          =   cell, a polysys cell containing a monomial basis for LT(Id)
%
% a             =   vector, contains indices of leading monomials
%
% n             =   scalar, number of variables
%
% CALLS
% -----
% 
% fite.m
%
% Kim Batselier, 2011-08-25, 2012-09-13: added 2nd syntax allowing for
% indices as input argument

if ~iscell(Asys)
    % the number of variables need to be given
    n=varargin{1};
    % convert indices to Asys element
    temp = cell(length(Asys),2);
    for i=1:length(Asys)
        temp{i,2} = fite(Asys(i),n);
    end
    clear Asys
    Asys = temp;
end

hasAll = 0;
n = size(Asys{1,2},2);
indices = zeros(1,n);

components = zeros(1,n);
counter = 1; % counts over Asys
indexcounter= 1; % counts over the indices

while ~isempty(find(components==0)) && counter <= size(Asys,1)
    
    if length(find(Asys{counter,2}))==1 && components(find(Asys{counter,2}))==0
        components(1,find(Asys{counter,2})) = Asys{counter,2}(find(Asys{counter,2})); 
        indices(indexcounter) = counter;
        indexcounter = indexcounter + 1;
    end
    counter = counter +1;        
    
end

hasAll = isempty(find(components==0));

end
