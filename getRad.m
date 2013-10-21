function radsys = getRad(polysys)
% radsys = getRad(polysys)
% ------------------------
% If a polynomial system polysys is not radical then calling this function
% will add polynomials such that the resulting polynomial system radsys is
% radical.
%
% radsys    =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations which are generators of the
%               radical ideal of polysys
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               original set of polynomial equations which are not
%               generators of a radical ideal.
%
% CALLS
% -----
% 
% Kim Batselier, 2011-01-23

warning off all

% number of unknowns
n = size(polysys{1,2},2);
% number of equations
n_eq = size(polysys,1);

radsys = polysys;

israd = 1;
for i = 1 : n

    % first need the monic generators
    [puni duni] = punivar(polysys,i);
    tol=eps*max(getMDim(polysys,duni))*20;
    puni=puni(1:find(abs(puni)>tol(end),1,'last')); % prune numerical zero higher order terms
    puni(abs(puni)<tol(end))=0; % set remaining numerical zeros to exact zero
    ppuni = prune(diffPoly(puni,1,1),1);
    ppuni = ppuni/norm(ppuni);
    [lcm h e] = getLCM(vec2polysys(puni,1),vec2polysys(ppuni,1),tol);
    h=h/norm(h);
    if  (length(h)~= length(puni)) || abs(puni*h')-1 > tol
        israd = 0;        
    end
    h=h(1:find(abs(h)>tol,1,'last'));
    h(abs(h)<tol)=0;
    
    temp = cell(1,2);
    coefI = find(h(1,:));
    temp{1,1} = h(coefI);    
    for j = 1 : length(coefI)
        temp{1,2} = [temp{1,2} ; [zeros(1,i-1) coefI(j)-1 zeros(1,n-i) ] ];
    end
    radsys{n_eq+i,1} = temp{1};
    radsys{n_eq+i,2} = temp{2};       
    
    
end

if israd
    clear radsys
    radsys = polysys;
end

end




