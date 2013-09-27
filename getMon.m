function fullBase = getMon(d,n,varargin)
% fullBase = getMon(d,n) or getMon(d,n,d0)
% ----------------------------------------------
%
% Returns a full canonical base of monomials of total degree d and in n
% variables. Each row of mon refers to a n-tuple of exponents of a monomial
% whereby each column corresponds to a variable.
%
% example: d= 2, n = 3
%
% fullBase =
% 
%      0     0     0
%      1     0     0
%      0     1     0
%      0     0     1
%      2     0     0
%      1     1     0
%      1     0     1
%      0     2     0
%      0     1     1
%      0     0     2
%
% CALLS
% -----
%
% getMonBase.m
%
% Kim Batselier, 2009-10
%
% Updated 2010-12-04: removed the monomial ordering argument. Replaced it
% by a minimal degree from which the monomial basis needs to be calculated.
% Further optimized the code by calling getMonBase once per for-iteration.

if ~isempty(varargin)
    d0 = varargin{1};
    if d0 > d
        fullBase = [];
        return
    end
else
    d0 = 0;
end

if d0 == 0
    fullBase = zeros(nchoosek(d+n,n),n);
    rowCounter = 2;
    for i = 1 : d,
        tempbase = getMonBase(i,n);
        fullBase(rowCounter:rowCounter+length(tempbase)-1,:)  = tempbase;
        rowCounter = rowCounter+length(tempbase);
    end
else
    rowCounter = 1;
    for i = d0 : d,
        tempbase = getMonBase(i,n);
        fullBase(rowCounter:rowCounter+length(tempbase)-1,:)  = tempbase;
        rowCounter = rowCounter+length(tempbase);
    end
end

end
