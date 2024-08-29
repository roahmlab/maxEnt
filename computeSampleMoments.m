function em = computeSampleMoments( data, l1hp )
%
%  data      -- k-by-dim array of data for which we wish to compute the
%               sample moments
%  l1hp      -- l-by-dim indexing array telling us where to look in the
%               moments matrix while building the hessian ( l < k )
% 
%  Computes the l sample moments from data 

em = zeros( size( l1hp, 1 ), 1 );
N = size( data, 1 );

em( 1 ) = 1;
for i = 2:size( l1hp, 1 )
    powers = repmat( l1hp( i, : ), N, 1 );
    em( i ) = 1/N * sum( prod( data.^powers, 2 ) );
end
