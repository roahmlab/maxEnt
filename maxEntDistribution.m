function [ lambda, entropy ] = maxEntDistribution( dim, num_moms, acc, em, ul, ll, l1hp, l12numhp, l1ap, l12numap, c_scale )
%
%  dim       -- distribution dimensionality
%  num_moms  -- nonnegative scalar integer
%  acc       -- nonnegative scalar integer for sparse gridding accuracy
%  em        -- l-by-1 vector of empirically computed moments for all
%                combinations with l1 norm less than num_moms
%  ul        -- 1-by-dim (the upper limit of integration)
%  ll        -- 1-by-dim (the lower limit of integration)
%  l1hp      -- l-by-dim indexing array telling us where to look in the
%               moments matrix while building the hessian
%  l12numhp  -- ( num_moms + 1 )^dim-by-1 indexing array telling us
%               where to look in the l1powers array while building the
%               hessian
%  l1ap      -- k-by-dim indexing array telling us where to look in the
%               moments matrix while building the hessian
%  l12numap  -- ( 2 * num_moms + 1 )^dim )-by-1 indexing array telling us
%               where to look in the l1powers array while building the
%               hessian
%  c_scale   -- 1-by-dim describing the scaling in each direction (useful
%               for plotting!)
%
% Computes the lambda_j for the maximum entropy distribution exp(\sum_{\abs{j} \leq n}
% \lambda_j x^j) such that its moments are equal to the emp_moments

%% parameters
% addpath( 'nwSpGr' );
numl1hp= size( l1hp, 1 );
numl1ap = size( l1ap, 1 );
numSteps = 100;
criterion = 1e-6;
mod = 1e-10;
plotting = 0; % only works for 3D data

%% preallocating nodes and powers
[ nodes, weights ] = nwspgr( 'KPU', dim, acc );
% need to apply change of variables on nodes
nodes = ( repmat( ul, size( nodes, 1 ), 1 ) - repmat( ll, size( nodes, 1 ), 1 ) ) .* nodes + repmat( ll, size( nodes, 1 ), 1 );
x = msspoly( 'x', dim );
[ ~, xn ] = isfree( x );
xhp = msspoly( [ numl1hp 1 ],  [ (1:numl1hp)' ones(numl1hp,1) ], repmat(xn',numl1hp,1), l1hp, ones(numl1hp,1) );
xap = msspoly( [ numl1ap 1 ],  [ (1:numl1ap)' ones(numl1ap,1) ], repmat(xn',numl1ap,1), l1ap, ones(numl1ap,1) );

%% newton's method for estimating lambda
% pre_lambda = randn( size( xhp, 1 ), 1 );
pre_lambda = zeros( size( xhp, 1 ), 1 );
moments = computeMoments( x, nodes', ul, ll, weights, xap, xhp, pre_lambda' );
[ H, am ] = moments2Hessian( num_moms, moments, l1hp, l12numap );
post_lambda = pre_lambda + pinv( H ) * ( em - am );
% sc = norm( post_lambda - pre_lambda, 1 )
sc = Inf;
counter = 1;

% foo = ones( size( em ) );
% foo( end ) = 1e5;
% em = foo .* em;
while ( sc > criterion && counter < numSteps );
    pre_lambda = post_lambda;
    moments = computeMoments( x, nodes', ul, ll, weights, xap, xhp, pre_lambda' );
    [ H, am ] = moments2Hessian( num_moms, moments, l1hp, l12numap );
    
    sc = norm( ( em - am )./em , Inf )
    
    if( any( isnan( H( : ) ) ) )
        break;
    end
    
    [ U, S, V ] = svd( H );
    if( cond( S ) > 1e5 || cond( S ) < 1e-5 )
        if( isinf( cond( S ) ) )
            break;
        end
        for i = 20:-10:-20
            M = cond( S ) * 0.8^i * eye( size( S, 1 ) );
            m2 = computeMoments( x, nodes', ul, ll, weights, xap, xhp, ( pre_lambda + cond( S ) * pinv( U * ( S + M ) * V' ) * ( em - am ) )' );
            [ ~, am2 ] = moments2Hessian( num_moms, m2, l1hp, l12numap );
            if( norm( em - am2, 2 ) < sc )
                %                 i
                break;
            end
        end
        post_lambda = pre_lambda + cond( S ) * 0.9^i * pinv( U * ( S + M ) * V' ) * ( em - am );
    else
        post_lambda = pre_lambda + pinv( U * S * V' ) * ( em - am );
    end
    
    %     sprintf( 'min diag = %f and cond = %f', min( abs( diag( S ) ) ), cond( S ) )
    %     h_mod = min( mod, 1e-3 * sc ); % as we get closer to the optimum we want the scaling to be less aggressive so that convergence happens quicker
    %     moddiag = h_mod * eye( size( S, 1  ) );
    %     moddiag( S > h_mod ) = 0;
    %     moddiag( : ) = 0;
    %     post_lambda = pre_lambda + pinv( U * ( S + moddiag ) * V' ) * ( em - am );
    
    counter = counter + 1;
end
if( any( isnan( H( : ) ) ) || isinf( cond( S ) ) )
    lambda = nan;
    entropy = nan;
else
    lambda = post_lambda;
    entropy = -sum( lambda .* am ) + log( abs( prod( c_scale ) ) );
end

% if( dim == 1 )
%     variance = am( 3 ) - 2 * am( 2 ) * am( 2 ) + am( 2 )^2 * am( 1 );
% end

%% plotting functionality (only works for 3D data at the moment)
% f = @(x) exp( lambda' * prod( repmat( x, size( l1hp, 1 ), 1 ).^l1hp, 2 ) );
% samplesfromDist = slicesample( ( ll + ul )/2 , 5000, 'pdf', f, 'thin', 5, 'burnin', 10000, 'width', 0.025 );
%
%
% figure( 1 );
% scatter3( samplesfromDist( :, 1 ) * c_scale( 1 ), samplesfromDist( :, 2 ) * c_scale( 2 ), ...
%     samplesfromDist( :, 3 ) * c_scale( 3 ) );

% % this plotting will work if the data is in 2D
% figure( 1 );
% cloudPlot( samplesfromDist( :, 1 ) * c_scale( 1 ), samplesfromDist( :, 2 ) * c_scale( 2 ) )
%
% % this ploting will work if the data is in 2D
% [ X, Y ] = meshgrid( ll(1):.01:ul(1), ll(2):0.01:ul(2) );
% exp_helper = exp( dmsubs( lambda' * xhp, x, [ X( : ) Y( : ) ]' ) );
% foo = reshape( exp_helper, size( X, 1 ), size( X, 2 ) );
% figure( 2 );
% mesh( X * c_scale( 1 ), Y * c_scale( 2 ), foo )


%
% if( plotting )
% %     [ X, Y, Z ] = sphere( 20 );
% %     exp_helper = exp( dmsubs( lambda' * xhp, x, [ X( : ) Y( : ) Z( : ) ]' ) );
% %     [ az, el, ~ ] = cart2sph( X( : ), Y( : ), Z( : ) );
% %     figure;
% %     PlotSphereIntensity( az * 180/pi, el * 180/pi, exp_helper' );
%      space = linspace( ll, ul, 1001 );
%      exp_helper = exp( dmsubs( lambda' * xhp, x, space ) );
%      figure; plot( space, exp_helper );
% end


