function [x,fval,exitflag] = spot_linprog(f,A,b,lb,ub)
    if nargin < 3
        error('Three arguments required.');
    end
    if nargin < 4, lb = []; end
    if nargin < 5, ub = []; end
    
    %[x,fval] = linprog(f,A,b,lb,ub);
    [x,fval] = glpk(f,A,b,lb,ub);
end