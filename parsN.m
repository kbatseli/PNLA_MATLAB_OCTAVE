function [N,tol]=parsN(polysys,d)

s=size(polysys,1);
n=size(polysys{1,2},2);

% normalize each coefficient vector to unit 1-norm such that the largest
% singular value of M(d) is upper bounded by sqrt(s)
for i=1:s
    polysys{i,1}=polysys{i,1}/norm(polysys{i,1},1);
    [Nf{i},nM,tol]=sparseN({polysys{i,1},polysys{i,2}},d);
end

temp=Nf{1};
for i=2:s
   inter=sparseN([temp Nf{i}],tol);
   temp = temp*inter(1:size(temp,2),:);
end

N=temp;


end