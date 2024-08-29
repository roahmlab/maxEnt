function moments = computeMoments( x, nodes, ul, ll, weights, xap, xhp, lambda )
%
%  x        -- dim-by-1 free msspoly
%  nodes    -- dim-by-k (arrangement of k precomputed nodes with limits taken into account in each dim)
%  ul       -- 1-by-dim (the upper limit of integration)
%  ll       -- 1-by-dim (the lower limit of integration)
%  weights  -- k-by-1 (weights of the k precomputed nodes)
%  xap       -- m-by-1 msspoly of all powers in x
%  xhp       -- n-by-1 msspoly of half powers in x
%  lambda   -- 1-by-n (list of all lambda)
% 
%  Computes all (instead of just half) moments of exp(\sum_{\abs{j} \leq num_moms} \lambda_j x^j)

moments = zeros( length( xap ), 1 );
exp_helper = exp( dmsubs( lambda * xhp, x, nodes ) );
det_cv = prod( ul - ll );

parfor i = 1:length( moments )
    mon_helper = dmsubs( xap( i ), x, nodes );
    moments( i ) = ( mon_helper .* exp_helper ) * weights * det_cv;
end
