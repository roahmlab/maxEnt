function q=recomp(x,p,M)
%  q=recomp(x,p,M)
%
%
%  Inverse of msspoly/decomp.
%
%  x -- n-by-1 free msspoly.
%  p -- m-by-n non-negative integers.
%  M -- N-by-m sparse matrix.
%
%  Constructs an N-by-1 msspoly p with entries:
%
%  p = M*repmat(x',m,1).^p;

if nargin < 2, error('2 inputs required.'); end


[free,xn] = isfree(x);
if ~free
    error('1st argument must be free.');
end
if size(x,2) ~= 1
    error('1st argument must be a column.');
end

if ~spot_isIntGE(p,-Inf)
    error('2nd arguments must be non-negative integers');
end

nx=size(x,1);
[mp,np]=size(p);

if nargin < 3, M = eye(mp); end


if size(p,1) ~= size(M,2)
    error('Powers and coefficients do not have matching dimensions.');
elseif length(x) ~= size(p,2)
    error('powers and variables do not have matching dimensions.');
end



N = size(M,1);
nz = reshape(find(M ~= 0),[],1);
[I,J] = ind2sub(size(M),nz);

vs = repmat(xn',length(J),1);


dim = [ N 1 ];
sub = [ I ones(size(I)) ];
var = vs;
pow = full(p(J,:));
coeff = reshape(full(M(nz)),[],1);

q = msspoly(dim,sub,var,pow,coeff);

end
