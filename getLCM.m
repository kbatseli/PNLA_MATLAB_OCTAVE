function [lcm h e] = getLCM(fsys,gsys,noiselevel)
% [lcm h e] = getLCM(fsys,gsys,tol)
% ---------------------------------
% Returns the approximate least common multiple of two given polynomials fsys and gsys. Also
% the factor h such that lcm = g*h is calculated.
%
% lcm       =   vector, coefficient vector of polynomial which is the least
%               common multiple of f and g
%
% h         =   vector, coefficient vector of polynomial which is the 
%               factor such that lcm = g*h
%
% e         =   scalar, mean squared error of least squares solution
%               lcm = g * h
%
% fsys      =   cell, polysys representation of the polynomial f
%
% gsys      =   cell, polysys representation of the polynomial g
%
% tol       =   scalar, tolerance: measure of how many digits are not
%               corrputed by noise. Default: eps
%
% CALLS
% -----
% 
% getM.m, getD0.m, getMex.m
%
% Kim Batselier, 2011-01-23, edited 2012-01-31: accepts polysys objects now
% warning off all
if nargin == 2
    noiselevel = 0;
end    

% check whether both f and g have same number of variables
nf = size(fsys{1,2},2);
ng = size(gsys{1,2},2);

if nf ~= ng
    error('Please provide 2 polynomials in the same number of variables');
end

% normalize coefficient vectors
fsys{1,1}=fsys{1,1}/norm(fsys{1,1});
gsys{1,1}=gsys{1,1}/norm(gsys{1,1});


% first get the degrees of the polynomials
df = getD0(fsys);
dg = getD0(gsys);

lcm = [];
h = [];
e = 0;
% need to start iteration over max(deg(f),deg(g))
d = max(df,dg);

% initialize the Macaulay matrices
Mf = getM(fsys,d,1);
[Q , ~, P] = qr(Mf','vector'); % very painful if c(d) is comparable to q(d)
Qf=Q(:,size(Mf,1));   % always of full row rank
Nf=Q(:,size(Mf,1)+1:end);

Mg = getM(gsys,d,1);
[Q , ~, P] = qr(Mg','vector'); % very painful if c(d) is comparable to q(d)
Ng = Q(:,size(Mg,1)+1:end); 
tol = 20*sum(size(Mf))*eps;


while isempty(lcm) && d <= df+dg
    
    if noiselevel % only for nonzero noiselevel we need to overwrite the tolerance
        tol=sqrt(2)*( cond(getM(fsys,d)) + cond(getM(gsys,d)) )*noiselevel;
    end
    [~, S Z] = svd(full(Ng'*Qf));
    
    if size(S,2)==1
        s=S(1);
    else
        diagS=diag(S);
        s=diagS(end);
    end
    
    if asin(s)<tol
        lcm = sparse(Qf*Z(:,end))';
        % decomposition of lcm in M via qr
        % Mg is always of full row rank!
        
        M=getM(gsys,d,1);        
        h =  (M'\lcm')';   
        e = norm(lcm-h*M);
    else
      d=d+1;
      [Qf Nf tol] = updateOrth(Qf,Nf,getMex(fsys,d,d-1,1),1);
      Ng = updateN(Ng,getMex(gsys,d,d-1,1),1);      
    end
end
    

if isempty(lcm)
    error('Method did not converge, try a smaller tolerance')
end

end
