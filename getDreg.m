function dreg = getDreg(polysys)
% dreg = getDreg(polysys)
% -----------------------
% Returns an upper bound for the degree at which the Hilbert function becomes the Hilbert
% polynomial of the given polynomial system
%
% dreg      =   scalar, regularity = degree for which M matrix has stable
%               nullity
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
%
% CALLS
% -----
%
% LITERATURE
% ----------
%
% Asymptotic Behaviour of the Degree of Regularity of Semi-Regular
% Polynomial Systems, Bardet M., J-C Faugere, B. Salvy, B-Y Yang
%
% Kim Batselier, 2010-04-21
dreg = 0;
for i = 1 : size(polysys,1)
    dreg = dreg + max(sum(polysys{i,2},2))-1;
end

dreg = dreg +1;

end
