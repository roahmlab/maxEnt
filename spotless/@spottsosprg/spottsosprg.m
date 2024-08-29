classdef spottsosprg < spotsqlprg
    properties
        sosExpr = {};
        sosTokens = {};
        tokenNo = 0;
        numSOS = 0;
    end
    methods (Access = private)
        function [flag,indet] = realPolyLinearInDec(pr,exp)
            [x,pow,Coeff] = decomp(exp);
            
            mtch = match(x,pr.variables);
            
            flag = ~(any(any(pow(:,mtch(mtch~=0)) > 1)) | ...
                     any(imag(Coeff(:)) ~= 0));
            
            if ~flag, 
                indet = [];
            else
                mtch = match(pr.variables,x);
                indet = x(find(mtch==0));
                if any(istrig(indet))
                    flag = 0;
                    indet = [];
                end
            end
        end
        
        function [flag,tIn,pIn] = trigPolyLinearInDec(pr,expr)
            [x,pow,Coeff] = decomp(expr);
            
            mtch = match(x,pr.variables);
            
            % TODO: conj. symmetric check.
            flag = ~(any(any(pow(:,mtch(mtch~=0)) > 1)));
            
            if ~flag, 
                indet = [];
            else
                mtch = match(pr.variables,x);
                indet = x(find(mtch==0));
                msk = istrig(indet);
                tIn = indet(find(msk));
                pIn = indet(find(~msk));
            end
        end
    end
    
    methods (Access = protected)
        function [pr,Q,phi,y,basis] = buildSOSDecomp(pr,expr)
            if ~spot_hasSize(expr,[1 1])
                error('buildSOSDecomp expects a scalar polynomial.');
            end

            decvar = pr.variables;

            [var,pow,M] = decomp(expr);
            mtch = match(var,decvar);
            b = 1:length(var);
            b(mtch(mtch ~= 0)) = [];
            indet = var(b);
            
            if length(indet) == 0
                [pr,Q] = pr.withPos(expr);
                phi = 1;
                return;
            end

            pow = pow(:,b);

            exponent_p_monoms = pow;
            csclasses={1:length(b)};
            exponent_m = monomialgeneration(exponent_p_monoms,csclasses);
    
            options = sdpsettings;
            options.verbose = 0;
            temp=sdpvar(1,1);
            tempops = options;
            tempops.solver = 'cdd,glpk,*';  % CDD is generally robust on these problems
            tempops.verbose = 0;
            tempops.saveduals = 0;
            [aux1,aux2,aux3,LPmodel] = export(set(temp>0),temp,tempops);  
            %disp('Reducing Monomials.');
            exponent_m = monomialreduction(exponent_m,exponent_p_monoms,options,csclasses,LPmodel);
            exponent_m = exponent_m{1};
    
            phi = recomp(indet,exponent_m,eye(size(exponent_m,1)));
            [pr,Q] = pr.newPSD(length(phi));
    
            decvar = [decvar ; mss_s2v(Q)];
            sosCnst = expr-phi'*Q*phi;

            %[pr,y,basis] = pr.withPolyEqs(expr-phi'*Q*phi);
            A = diff(sosCnst,decvar);
            b = subs(sosCnst,decvar,0*decvar);
            [var,pow,Coeff] = decomp([b A].');
            
            [pr,y] = pr.withEqs(Coeff'*[1;decvar]);
            basis = recomp(var,pow,eye(size(pow,1)));
        end
    end
    
    methods
        function pr = spottsosprg(varargin)
            pr@spotsqlprg(varargin{:});
        end
        
        function [pr,tokens] = withSOS(pr,expr,c,s,mults)
            if ~pr.realPolyLinearInDec(expr)
                error(['Coefficients must be real, indeterminates ' ...
                       'non-trigonometric, and expression must ' ...
                      'be linear in decision variables.']);
            end
            tokens = pr.tokenNo + reshape(1:prod(size(expr)),size(expr));
            pr.tokenNo = max(tokens);
            pr.numSOS = pr.numSOS + prod(size(expr));
            
            if nargin == 2
                pr.sosExpr{end+1} = {expr};
            elseif nargin == 5
                pr.sosExpr{end+1} = {expr,c,s,mults};
            else
                error('Need 2 or 5 arguments.');
            end
            pr.sosTokens{end+1} = tokens;
        end
        
        function [pr,y,basis] = withPolyEqs(pr,expr)
            if ~pr.realPolyLinearInDec(expr)
                error(['Coefficients must be real, indeterminates ' ...
                       'non-trigonometric, and expression must ' ...
                       'be linear in decision variables.']);
            end

            expr = expr(:);
            decvar = pr.variables;
            
            [indet,pow,M] = decomp(expr,decvar);
            
            monom = recomp(indet,pow,eye(size(pow,1)));
            
            [I,J,S] = find(M);
            
            [pr,y] = pr.withEqs(S);
            
            basis = monom(J);
        end
        
        function [pr,tokens] = withSOSMatrix(pr,expr)
            [lindec,indet] = pr.realPolyLinearInDec(expr);
            if ~lindec
                error(['Coefficients must be real, indeterminates ' ...
                       'non-trigonometric, and expression must ' ...
                      'be linear in decision variables.']);
            end
            
            if size(expr,1) ~= size(expr,2) || size(expr,1) == 0
                error('Expression must be a square non-empty matrix.');
            end
            
            x = msspoly('x',size(expr,1));
            expr = anonymize(expr,'y',indet);
            [pr,tokens] = withSOS(pr,x'*expr*x);
        end
        
        function [pr,tokens] = withUniTrigSOSMatrix(pr,expr)
            if size(expr,1) ~= size(expr,2) || size(expr,1) == 0
                error('Expression must be a square non-empty matrix.');
            end
            
            [lindec,trigIndet,polyIndet] = pr.trigPolyLinearInDec(expr);
            if ~lindec
                error(['Expression ' ...
                       'must be linear in decision variables.']);
            elseif ~isempty(polyIndet) || (length(trigIndet) ~= 1)
                error('Only a single, trigonometric indeterminate allowed.');
            end
            
            
            % How to handle these?  Special case? Maybe...

        end
        
        function [pr,poly,coeff] = newFreePoly(pr,basis,n)
            if nargin < 3, n = 1; end
            [pr,coeff] = pr.newFree(length(basis)*n);
            poly = reshape(coeff,n,length(basis))*basis;
        end
        
        function sol = optimize(pr,varargin)
            if nargin < 2, objective = msspoly(0); end

            Q = cell(pr.numSOS,1);
            phi = cell(pr.numSOS,1);
            y   = cell(pr.numSOS,1);
            basis   = cell(pr.numSOS,1);
            
            i = 1;
            for k = 1:length(pr.sosExpr)
                expr = pr.sosExpr{k};
                
                if spot_hasSize(expr,[1 1])
                    expr = expr{1};
                    for j = (1:prod(size(expr)))
                        [pr,Q{i},phi{i},y{i},basis{i}] = pr.buildSOSDecomp(expr(j));
                        i = i+1;
                    end
                else
                    c = expr{2}; s = expr{3}; mults = expr{4};
                    expr = expr{1};
                    
                    for j = 1:prod(size(expr))
                        [pr,Q{i},phi{i},y{i},basis{i}] = pr.buildTrigSOSDecomp(expr(j),c,s,mults);
                        i = i+1;
                    end
                end
            end

            sqlsol = optimize@spotsqlprg(pr,varargin{:});

            
            sol = spotsossol(sqlsol,Q,phi,y,basis);
        end
    end
end