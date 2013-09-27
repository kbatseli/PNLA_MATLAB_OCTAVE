function tau = tauBound(epsilon,psys,qsys)
% tau = tauBound(epsilon,psys,qsys)
% ---------------------------------
% 
% Computes an estimate for tolerance tau to use to determine whether principal angle is 0.
%
% tau         = scalar, tolerance 
%
% epsilon     = scalar, upper bound for noise level
%
% psys        = polysys object of multivariate polynomial
%
% qsys        = polysys object of multivariate polynomial
%
% REFS
% ----
%
% A geometrical approach to finding multivariate approximate LCMs and GCDs
%
% % Kim Batselier, 2012-07-12

n = size(psys{1,2},2);
m = zeros(1,2); %contains constant coefficients of psys and qsys

constantexponent = zeros(1,n);

posfound = 0;
teller = 1;
while ~posfound & teller <= size(psys{1,2},1)
	if sum(psys{1,2}(teller,:)==constantexponent)==n
		posfound = 1;
		m(1) = psys{1,1}(teller);
	else
		teller = teller + 1;
	end
end


posfound = 0;
teller = 1;
while ~posfound & teller<=size(qsys{1,2},1)
	if sum(qsys{1,2}(teller,:)==constantexponent)==n
		m(2) = qsys{1,1}(teller);
		posfound = 1;
	else
		teller = teller + 1;
	end
end	

% check for nontrivial case
if 2*abs(m(1)) < norm(psys{1,1},1) || 2*abs(m(2)) < norm(qsys{1,1},1)
	disp('Trivial case: bound does not apply');
	return
end

tau = epsilon*(1/(2*abs(m(1)) -norm(psys{1,1},1)) + 1/(2*abs(m(2)) - norm(qsys{1,1},1)) );
end
