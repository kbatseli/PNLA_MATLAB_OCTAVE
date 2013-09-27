function polysys = lti2polysys(strings,y,u,varargin)
% polysys = lti2polysys(strings,y,u,'na',na,'nb',nb,...)
% ------------------------------------------------------
% 
% Converts a cell of strings, representing a polynomial system from a PEM
% of a LTI model, into a polysys cell.
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% strings   =   cell, each cell contains 1 string from a LTI model
%
% y         =   vector, contains the measured output for system
%               identification
%
% u         =   vector, contains the measured input for system
%               identification
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

for k = 1 : 2 : length(varargin)
    if varargin{k} == 'na'
        na = varargin{k+1};
    elseif varargin{k} == 'nb'
        nb = varargin{k+1};
    elseif varargin{k} == 'nc'
        nc = varargin{k+1};
    elseif varargin{k} == 'nc'
        nd = varargin{k+1};        
    elseif varargin{k} == 'nf'
        nf = varargin{k+1};
    elseif varargin{k} == 'nl'
        nl = varargin{k+1}; 
    elseif varargin{k} == 'ne'
        ne = varargin{k+1};         
    end
end

n = na+nb+nc+nd+nf+nl+ne; % number of unknowns

polysys = cell(length(strings),2);

for stringcounter = 1 : length(strings) % process each string
    
    temp = strings{stringcounter};
    l = length(temp);
    
    pindex = strfind(temp,'+'); % find the +'s
    mindex = strfind(temp,'-'); % find the -'s
    
    indices = sort([pindex mindex]); %borders of the terms  
    
   
    if isempty(indices)
        indices = [1];        
    elseif indices(1) ~= 1
        % first term does not start with '-'
        indices = [1 indices];
    end
    
    indices = [indices l+1];
        
    coefs = zeros(1,length(indices)-1);
    polysys{stringcounter,2} = [];
    for j = 1 : length(indices)-1   % for each term
        term = temp(indices(j):indices(j+1)-1); % isolate the term
                
        % process the coefficients        
        prodindex = [0 strfind(term,'*') length(term)+1];
        tempcoef = '1';
        exp = zeros(1,n);
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
            
%             if factor(1) == '-'
                switch factor(checkindex)
                    case 'u'                        
                        blindex = strfind(factor,'[');
                        brindex = strfind(factor,']');
                        expind = strfind(factor,'^');
                        expvalue = '1';
                        if ~isempty(expind)
                            expvalue = factor(expind+1:end);                            
                        end
                        tempcoef = [tempcoef '* u(' factor(blindex+1:brindex-1) ')^' expvalue];
                    case 'y'
                        blindex = strfind(factor,'[');
                        brindex = strfind(factor,']');
                        expind = strfind(factor,'^');
                        expvalue = '1';
                        if ~isempty(expind)
                            expvalue = factor(expind+1:end);                            
                        end
                        tempcoef = [tempcoef '* y(' factor(blindex+1:brindex-1) ')^' expvalue];
                    case 'a'
                        %find out which a
                        blindex = strfind(factor,'[');
                        brindex = strfind(factor,']');
                        aindex = str2num(factor(blindex+1:brindex-1));
                        % check for exponent
                        expind = strfind(factor,'^');
                        expvalue = 1;
                        if ~isempty(expind)
                            expvalue = str2num(factor(expind+1:end));                            
                        end
                        exp = exp + [zeros(1,aindex-1) expvalue zeros(1,n-aindex)];
                    case 'b'
                        blindex = strfind(factor,'[');
                        brindex = strfind(factor,']');
                        bindex = str2num(factor(blindex+1:brindex-1));
                        % check for exponent
                        expind = strfind(factor,'^');
                        expvalue = 1;
                        if ~isempty(expind)
                            expvalue = str2num(factor(expind+1:end));
                        end
                        exp = exp + [zeros(1,na) zeros(1,bindex-1) expvalue zeros(1,n-na-bindex)];
                    case 'c'
                        blindex = strfind(factor,'[');
                        brindex = strfind(factor,']');
                        cindex = str2num(factor(blindex+1:brindex-1));
                        % check for exponent
                        expind = strfind(factor,'^');
                        expvalue = 1;
                        if ~isempty(expind)
                            expvalue = str2num(factor(expind+1:end));
                        end                        
                        exp = exp + [zeros(1,na+nb) zeros(1,cindex-1) expvalue zeros(1,n-na-nb-cindex)];
                    case 'd'
                        blindex = strfind(factor,'[');
                        brindex = strfind(factor,']');
                        dindex = str2num(factor(blindex+1:brindex-1));
                        % check for exponent
                        expind = strfind(factor,'^');
                        expvalue = 1;
                        if ~isempty(expind) % exponent bigger than 1
                            expvalue = str2num(factor(expind+1:end)); 
                        end                        
                        exp = exp + [zeros(1,na+nb+nc) zeros(1,dindex-1) expvalue zeros(1,n-na-nb-nc-dindex)];
                    case 'f'
                        blindex = strfind(factor,'[');
                        brindex = strfind(factor,']');
                        findex = str2num(factor(blindex+1:brindex-1));
                        % check for exponent
                        expind = strfind(factor,'^');
                        expvalue = 1;
                        if ~isempty(expind) % exponent bigger than 1
                            expvalue = str2num(factor(expind+1:end)); 
                        end                        
                        exp = exp + [zeros(1,na+nb+nc+nd) zeros(1,findex-1) expvalue zeros(1,n-na-nb-nc-nd-findex)];
                    case 'l'
                        blindex = strfind(factor,'[');
                        brindex = strfind(factor,']');
                        lindex = str2num(factor(blindex+1:brindex-1));
                        % check for exponent
                        expind = strfind(factor,'^');
                        expvalue = 1;
                        if ~isempty(expind) % exponent bigger than 1
                            expvalue = str2num(factor(expind+1)); 
                        end                        
                        exp = exp + [zeros(1,na+nb+nc+nd+nf) zeros(1,lindex-1) expvalue zeros(1,n-na-nb-nc-nd-nf-lindex)];
                    case 'e'
                        blindex = strfind(factor,'[');
                        brindex = strfind(factor,']');
                        eindex = str2num(factor(blindex+1:brindex-1));
                        % check for exponent
                        expind = strfind(factor,'^');
                        expvalue = 1;
                        if ~isempty(expind) % exponent bigger than 1
                            expvalue = str2num(factor(expind+1)); 
                        end                        
                        exp = exp + [zeros(1,na+nb+nc+nd+nf+nl) zeros(1,eindex-1) expvalue zeros(1,n-na-nb-nc-nd-nf-nl-eindex)];                          
                    otherwise
                        tempcoef = [tempcoef '* ' '(' factor(checkindex:end) ')' ];
                end

            
        end
        
        coefs(j)= eval(tempcoef);
        polysys{stringcounter,2} = [polysys{stringcounter,2};exp];
    end
    polysys{stringcounter,1} = coefs;
    
    
    

end
