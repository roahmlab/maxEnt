function [ H, am ] = moments2Hessian( num_moms, moments, l1hp, l12numap )
%
%  num_moms  -- nonnegative scalar integer
%  moments   -- k-by-1 nonnegative array with each row describing all
%               powers for each dimension remember moments must contain all
%               possible combination of powers up to 2*num_moms
%  l1hp      -- l-by-dim indexing array telling us where to look in the
%               moments matrix while building the hessian ( l < k )
%  l12numap  -- ( 2 * num_moms + 1 )^dim )-by-1 indexing array telling us
%               where to look in the l1powers array while building the 
%               hessian
% 
%  Computes the Hessian (H) from moments vector and the analytical moments (am) up
%  to num_moms + 1

numl1hp = size( l1hp, 1 );
H = zeros( numl1hp );
am = zeros( numl1hp, 1 );

for i = 1:numl1hp
    for j = 1:numl1hp
        helper = l12numap( base2decnonstring( l1hp( i, : ) + l1hp( j, : ), 2 * num_moms + 1 ) + 1 );
        H( i, j ) = moments( helper + 1 );
    end
    am( i ) = moments( l12numap( base2decnonstring( l1hp( i, : ), 2 * num_moms + 1 ) + 1 ) + 1 );
end
