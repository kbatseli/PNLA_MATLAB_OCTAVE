function matrix2tex(M,n)
% matrix2tex(M,n)
% ---------------
%
% Converts a Matlab matrix to a latex file using bordermatrix with all
% columns labeled.
%
% CALLS
% -----
%
% Kim Batselier, 2012-05-14

[p q] = size(M);

d = sum(frte(n,q));


openingstring = ['\bordermatrix{\text{}'];
monomial = '&1';

for i = 2:q
    exponent = frte(n,i);
    monomial = [monomial '&'];
    for j = 1:n
        if exponent(j)
            if exponent(j)==1
                monomial = [monomial 'x_' num2str(j)]; 
            else
                monomial = [monomial 'x_' num2str(j) '^' num2str(exponent(j))];
            end
        end
    end    
end

openingstring = [openingstring monomial '\cr'];

 fid = fopen('matrix.txt','w');
 fprintf(fid,'%s\n',openingstring);
 
 numberstring =[];
 for i = 1 : p
     numberstring = [];
     for j = 1:q
         numberstring = [numberstring '&' num2str(M(i,j))];
     end
     numberstring = [numberstring '\cr'];
     fprintf(fid,'%s\n',numberstring);
 end
fprintf(fid,'%s\n\r','}');
 fclose(fid);


end