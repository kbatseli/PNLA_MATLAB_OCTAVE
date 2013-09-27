function V = evalV(Vstring,y,u,roots,varargin)
%  V = evalV(Vstring,y,u,roots,varargin)
% --------------------------------------
% 
% Evaluates the given cost function, given as the string Vstring, in the
% roots of the polynomial system.
%
% V         =   vector, contains the evaluated cost function
%
% y         =   vector, contains the measured output for system
%               identification
%
% u         =   vector, contains the measured input for system
%               identification
%
% roots     =   matrix, each column corresponds with a solution of the
%               polynomial system
%
% na        =   scalar, optional: number of 'a' parameters, other optional
%               parameters are nb, nc, nd, nf, nl (lagrange multipliers)
%               and ne (prediction errors).
%
% CALLS
% -----
%
% Kim Batselier, 2011-10

na=0;
nb=0;
nc=0;
nd=0;
nf=0;
nl=0;
ne=0;

for i = 1 : 2 : length(varargin)
    if varargin{i} == 'na'
        na = varargin{i+1};
    elseif varargin{i} == 'nb'
        nb = varargin{i+1};
    elseif varargin{i} == 'nc'
        nc = varargin{i+1};
    elseif varargin{i} == 'nc'
        nd = varargin{i+1};        
    elseif varargin{i} == 'nf'
        nf = varargin{i+1};
    elseif varargin{i} == 'nl'
        nl = varargin{i+1};
    elseif varargin{i} == 'ne'
        ne = varargin{i+1};           
    end
end

% roots zijn rijen van eerste graad uit kernel
[n nRoots] = size(roots);

V = zeros(1,nRoots);

l = length(Vstring);

pindex = strfind(Vstring,'+'); % find the +'s
mindex = strfind(Vstring,'-'); % find the -'s

indices = sort([pindex mindex]); %borders of the terms    

if indices(1) ~= 1
    % first term does not start with '-'
    indices = [1 indices];
end

indices = [indices l+1];

for j = 1 : length(indices)-1   % for each term
    
    term = Vstring(indices(j):indices(j+1)-1);
    prodindex = [0 strfind(term,'*') length(term)+1];    
    tempcoef = '1';
    
    for k = 1 : length(prodindex)-1   % each factor
        
        factor = term(prodindex(k)+1:prodindex(k+1)-1);
        checkindex = 1; % default index for k ~= 1
        signchar = '';
        if factor(1) == '-'                
            checkindex = 2; % first char is a sign, start to check from second
            tempcoef = [tempcoef '*(-1)'];
        elseif factor(1) == '+'
            checkindex = 2; % first char is a sign, start to check from second
        end
        
        switch factor(checkindex)
            case 'u'
                blindex = strfind(factor,'[');
                brindex = strfind(factor,']');
                % check for exponent
                expind = strfind(factor,'^');
                expvalue = 1;
                if ~isempty(expind)
                    expvalue = str2num(factor(expind+1:end));                            
                end                
                tempcoef = [tempcoef '.* u(' factor(blindex+1:brindex-1) ')^' num2str(expvalue)];
            case 'y'
                blindex = strfind(factor,'[');
                brindex = strfind(factor,']');
                expind = strfind(factor,'^');
                expvalue = 1;
                if ~isempty(expind)
                    expvalue = str2num(factor(expind+1:end));                            
                end
                tempcoef = [tempcoef '.* y(' factor(blindex+1:brindex-1) ')^' num2str(expvalue)];
            case 'a'
                %find out which a
                blindex = strfind(factor,'[');
                brindex = strfind(factor,']');
                aindex = str2num(factor(blindex+1:brindex-1));
                rowindex = aindex;
                % check for exponent
                expind = strfind(factor,'^');
                expvalue = 1;
                if ~isempty(expind)
                    expvalue = str2num(factor(expind+1:end));                            
                end
                tempcoef = [tempcoef '.* roots(' num2str(rowindex)  ',:).^' num2str(expvalue)];
            case 'b'
                blindex = strfind(factor,'[');
                brindex = strfind(factor,']');
                bindex = str2num(factor(blindex+1:brindex-1));
                rowindex = na+bindex;
                % check for exponent
                expind = strfind(factor,'^');
                expvalue = 1;
                if ~isempty(expind)
                    expvalue = str2num(factor(expind+1:end));
                end
                tempcoef = [tempcoef '.* roots(' num2str(rowindex)  ',:).^' num2str(expvalue)];
            case 'c'
                blindex = strfind(factor,'[');
                brindex = strfind(factor,']');
                cindex = str2num(factor(blindex+1:brindex-1));
                rowindex = na+nb+cindex;
                % check for exponent
                expind = strfind(factor,'^');
                expvalue = 1;
                if ~isempty(expind)
                    expvalue = str2num(factor(expind+1:end));
                end
                tempcoef = [tempcoef '.* roots(' num2str(rowindex)  ',:).^' num2str(expvalue)];
            case 'd'
                blindex = strfind(factor,'[');
                brindex = strfind(factor,']');
                dindex = str2num(factor(blindex+1:brindex-1));
                rowindex = na+nb+nc+dindex;
                % check for exponent
                expind = strfind(factor,'^');
                expvalue = 1;
                if ~isempty(expind) % exponent bigger than 1
                    expvalue = str2num(factor(expind+1:end)); 
                end
                tempcoef = [tempcoef '.* roots(' num2str(rowindex)  ',:).^' num2str(expvalue)];
            case 'f'
                blindex = strfind(factor,'[');
                brindex = strfind(factor,']');
                findex = str2num(factor(blindex+1:brindex-1));
                rowindex = na+nb+nc+nd+findex;
                % check for exponent
                expind = strfind(factor,'^');
                expvalue = 1;
                if ~isempty(expind) % exponent bigger than 1
                    expvalue = str2num(factor(expind+1:end)); 
                end
                tempcoef = [tempcoef '.* roots(' num2str(rowindex)  ',:).^' num2str(expvalue)];
            case 'l'
                blindex = strfind(factor,'[');
                brindex = strfind(factor,']');
                lindex = str2num(factor(blindex+1:brindex-1));
                rowindex = na+nb+nc+nd+nf+lindex;
                % check for exponent
                expind = strfind(factor,'^');
                expvalue = 1;
                if ~isempty(expind) % exponent bigger than 1
                    expvalue = str2num(factor(expind+1)); 
                end
                tempcoef = [tempcoef '.* roots(' num2str(rowindex)  ',:).^' num2str(expvalue)];    
            case 'e'
                blindex = strfind(factor,'[');
                brindex = strfind(factor,']');
                eindex = str2num(factor(blindex+1:brindex-1));
                rowindex = na+nb+nc+nd+nf+nl+eindex;
                % check for exponent
                expind = strfind(factor,'^');
                expvalue = 1;
                if ~isempty(expind) % exponent bigger than 1
                    expvalue = str2num(factor(expind+1)); 
                end
                tempcoef = [tempcoef '.* roots(' num2str(rowindex)  ',:).^' num2str(expvalue)];                 
            otherwise
                tempcoef = [tempcoef '.* ' '(' factor(checkindex:end) ')' ];
        end            
    end
    
    % update V
    V = V + eval(tempcoef);
end

end