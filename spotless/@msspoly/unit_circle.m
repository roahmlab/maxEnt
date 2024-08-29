function [c,s] = unit_circle(z)
%
% [c,s] = unit_circle(z)
%
% z -- n-by-1 free msspoly
%
% c -- (z+z.^-1)/2
% s -- (z-z.^-1)/(2j)
    
    n = size(z,1);
    [f,xn] = isfree(z);
    
    if ~f || size(z,2) ~= 1 
        error('First argument must be n-by-1 free msspoly.');
    end

    c = msspoly(size(z),...
                [[(1:n) (1:n)]' ones(2*n,1)],...
                [ xn ; xn ],...
                [ ones(n,1) ; -ones(n,1) ],...
                0.5*ones(2*n,1));
    s = msspoly(size(z),...
                [[(1:n) (1:n)]' ones(2*n,1)],...
                [ xn ; xn ],...
                [ ones(n,1) ; -ones(n,1) ],...
                0.5*j*[-ones(n,1);ones(n,1)]);
end