function [lcm h e] = getLCM(fsys,gsys,noiselevel,notsparse)
% [lcm h e] = getLCM(fsys,gsys,noiselevel,notsparse)
% --------------------------------------------------
% Returns the approximate least common multiple of two given polynomials fsys and gsys. Also
% the factor h such that lcm = g*h is calculated.
%
% lcm       	=   vector, coefficient vector of polynomial which is the least
%               	common multiple of f and g
%
% h         	=   vector, coefficient vector of polynomial which is the 
%               	factor such that lcm = g*h
%
% e         	=   scalar, mean squared error of least squares solution
%               	lcm = g * h
%
% fsys      	=   cell, polysys representation of the polynomial f
%
% gsys      	=   cell, polysys representation of the polynomial g
%
% noiselevel	=   scalar, upper bound on the error on the coefficients
%					of fsys and gsys. Default: tolerance for numerical rank Macaulay matrix
%
% notsparse     =   boolean, if set to 0, then a sparse rank revealing QR
%                   is used. If set to 1, a dense SVD is used instead.
%                   Default is 0.
%
% CALLS
% -----
% 
% getM.m, getD0.m, getMex.m, updateOrth.m, updateN.m
%
% Kim Batselier, 2011-01-23, edited 2012-01-31: accepts polysys objects now
% warning off all
if nargin == 2
    noiselevel = 0;
    notsparse = 0;
elseif nargin == 3
	notsparse = 0;
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

% initialization Macaulay matrix and orthogonal basis kernel
if ~notsparse
    Mf = getM(fsys,d,1);
    [Q , R, P] = qr(Mf','vector'); % very painful if c(d) is comparable to q(d)
    Qf=Q(:,size(Mf,1));   % always of full row rank
	Nf=Q(:,size(Mf,1)+1:end);

	Mg = getM(gsys,d,1);
	[Q ,R, P] = qr(Mg','vector'); % very painful if c(d) is comparable to q(d)
	Ng = Q(:,size(Mg,1)+1:end); 
	tol = 20*sum(size(Mf))*eps;
else
    Mf=getM(fsys,d)';
    [U S V]=svd(Mf);
    s=diag(S);
    tol=max(size(Mf))*eps(s(1));
    r=sum(s > tol );
    Qf=U(:,1:r);
    Nf=U(:,r+1:end);
    
    Mg = getM(gsys,d);
    Ng=null(Mg);
end


while isempty(lcm) && d <= df+dg
    
    if noiselevel % only for nonzero noiselevel we need to overwrite the tolerance
        tol=sqrt(2)*( cond(getM(fsys,d)) + cond(getM(gsys,d)) )*noiselevel;
    end
    [Y, S Z] = svd(full(Ng'*Qf));
    
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
        [asin(s) d tol]
        M=getM(gsys,d,1);        
        h =  (M'\lcm')';   
        e = norm(lcm-h*M);
    else
      d=d+1;
      if ~notsparse
		[Qf Nf tol] = updateOrth(Qf,Nf,getMex(fsys,d,d-1,1),1);
		Ng = updateN(Ng,getMex(gsys,d,d-1,1),1);        
	  else
		[Qf Nf tol] = updateOrth(Qf,Nf,getMex(fsys,d,d-1,1));
		Ng = updateN(Ng,getMex(gsys,d,d-1,1));      
	  end
      
    end
end
    

if isempty(lcm)
    error('Method did not converge, try a smaller tolerance')
end

end
