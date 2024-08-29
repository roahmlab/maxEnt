function d = base2decnonstring(h,b)
%BASE2DEC Convert base B string to decimal integer.
%   BASE2DEC(S,B) converts the string number S of base B into its 
%   decimal (base 10) equivalent.  B must be an integer
%   between 2 and 36. S must represent a non-negative integer value.
%
%   If S is a character array, each row is interpreted as a base B string.
%
%   Example
%      base2dec('212',3) returns 23
%
%   See also DEC2BASE, HEX2DEC, BIN2DEC.

%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 1.17.4.9 $  $Date: 2011/02/15 00:54:07 $

%   Douglas M. Schwarz, 18 February 1996

[m,n] = size(h);
bArr = [ones(m,1) cumprod(b(ones(m,n-1)),2)];
d = sum( ( bArr .* fliplr( h ) ), 2 );
