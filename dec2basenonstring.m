function s = dec2basenonstring(d,b,nin)
%DEC2BASENONSTRING Convert decimal integer to base B.
%   DEC2BASENONSTRING(D,B) returns the representation of D as a string in
%   base B.  D must be a non-negative integer array smaller than 2^52
%   and B must be an integer.
%
%   DEC2BASENONSTRING(D,B,N) produces a representation with at least N digits.
%
%   Examples
%       dec2base(23,3) returns 212
%       dec2base(23,3,5) returns 00212

% Original by Douglas M. Schwarz, Eastman Kodak Company, 1996.
if nargin<2
    error(nargchk(2,3,nargin,'struct'));
end
d = d(:);
if ~(isnumeric(d) || ischar(d)) || any(d ~= floor(d)) || any(d < 0) || any(d > 1/eps)
    error(message('MATLAB:dec2basenonstring:FirstArg'));
end
if ~isscalar(b) || ~(isnumeric(b) || ischar(b)) || b ~= floor(b)
    error(message('MATLAB:dec2basenonstring:SecondArg'));
end
if nargin == 3
    if ~isscalar(nin) || ~(isnumeric(nin) || ischar(nin)) || nin ~= floor(nin) || nin < 0
        error(message('MATLAB:dec2basenonstring:ThirdArg'));
    end
end
d = double(d);
b = double(b);
n = max(1,round(log2(max(d)+1)/log2(b)));
while any(b.^n <= d)
    n = n + 1;
end
if nargin == 3
    n = max(n,nin);
end
s(:,n) = rem(d,b);
% any(d) must come first as it short circuits for empties
while any(d) && n >1
    n = n - 1;
    d = floor(d/b);
    s(:,n) = rem(d,b);
end
