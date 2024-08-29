%% load data
addpath( '/nwSpGr' );
% choose the scale so that the data really fills up the range
percentile = 97;

%% process and scale data
NSamples = 10000;
predata = rand( NSamples, 3 );
predata( :, 1 ) = predata( :, 1 ) .* ( 5 - ( -5 ) ) + ( -5 );
predata( :, 2 ) = predata( :, 2 ) .* ( 0.33 - ( -0.33 ) ) + ( -0.33 );
predata( :, 3 ) = predata( :, 3 ) .* ( 0.001 - ( -0.001 ) ) + ( -0.001 );
[ data, postdata, coord_scale ] = scaleDataPER( predata, percentile );

%% parameters
dim = 3;
num_moms = 4; % need to compute twice as many of these in order to fill the moment matrix!
accuracy = 8; % has to do with numerical accuracy of moment computation ( generally should be larger than the number of moments + 1 )
% integration limits
ll = [ -1 -1 -1 ];
ul = [ 1 1 1 ];
idx = all( [ data( :, 1 ) <= ul( 1 ) data( :, 2 ) <= ul( 2 ) data( :, 3 ) <= ul( 3 ) ], 2 ) & ...
    all( [ data( :, 1 ) >= ll( 1 ) data( :, 2 ) >= ll( 2 ) data( :, 3 ) >= ll( 3 ) ], 2 );
data = data( idx, : );
predata = predata( idx, : );
pll = min( predata );
pul = max( predata );
[ sphdata( :, 1 ), sphdata( :, 2 ), sphdata( :, 3 ) ] = cart2sph( predata( :, 1 ), predata( :, 2 ), predata( :, 3 ) );
polll = min( sphdata );
polul = max( sphdata );

%% precompute moment indices to go in either direction (NEED THIS INFORMATION TO COMPUTE EMPIRICAL MOMENTS IN RIGHT ORDER)
hp = dec2basenonstring( ( 1:( ( num_moms + 1 )^dim ) ) - 1, num_moms + 1, dim );
num2l1hp = find( sum( hp, 2 ) <= num_moms );
l1hp = hp( num2l1hp, : );
numl1hp= size( l1hp, 1 );
l12numhp = -1 * ones( length( hp ), 1 );
l12numhp( num2l1hp ) = 0:( numl1hp - 1 );

% transform whole array of numbers into base 2 * num_moms + 1 with dim digits
ap = dec2basenonstring( ( 1:( ( 2 * num_moms + 1 )^dim ) ) - 1, 2 * num_moms + 1, dim );
num2l1ap = find( sum( ap, 2 ) <= 2 * num_moms );
l1ap = ap( num2l1ap, : );
numl1ap = size( l1ap, 1 );
l12numap = -1 * ones( length( ap ), 1 );
l12numap( num2l1ap ) = 0:( numl1ap - 1 );

%% go from data to empirical moments (em)
em = computeSampleMoments( data, l1hp );

%% use Newtons method to go from empirical moments to a distribution that satisfies requirements
[ lambda, entropy ] = maxEntDistribution( dim, num_moms, accuracy, em, ul, ll, l1hp, l12numhp, l1ap, l12numap, coord_scale );

%% rescale the lambdas
% need to plot p(r) going from p(x,y,z)
x = msspoly( 'x', dim );
[ ~, xn ] = isfree( x );
xhp = msspoly( [ numl1hp 1 ],  [ (1:numl1hp)' ones(numl1hp,1) ], repmat(xn',numl1hp,1), l1hp, ones(numl1hp,1) );
xap = msspoly( [ numl1ap 1 ],  [ (1:numl1ap)' ones(numl1ap,1) ], repmat(xn',numl1ap,1), l1ap, ones(numl1ap,1) );
lambdasc = lambda .* dmsubs( xhp, x, 1./coord_scale' ); % the lambdas for the scaled distribution

numl1hp= size( l1hp, 1 );
numl1ap = size( l1ap, 1 );
[ nodes, weights ] = nwspgr( 'KPU', dim - 1, 25 );

%% creating the distribution for original data
pxlims = pll( 1 ):0.01:pul( 1 );
prex = zeros( size( pxlims ) );
sl = pll( 2:3 );
su = pul( 2:3 );
det_s = prod( su - sl );
spnodes = ( repmat( su, size( nodes, 1 ), 1 ) - repmat( sl, size( nodes, 1 ), 1 ) ) .* nodes + repmat( sl, size( nodes, 1 ), 1 );
parfor i = 1:length( prex )
    eval_sp = zeros( size( spnodes, 1 ), 3 );
    eval_sp( :, 1 ) = pxlims( i );
    eval_sp( :, 2 ) = spnodes( :, 1 );
    eval_sp( :, 3 ) = spnodes( :, 2 );
    prex( i ) = sum( exp( dmsubs( lambdasc' * xhp, x, eval_sp' ) )' .* weights .* det_s ) * 1/prod( coord_scale );
end

pylims = pll( 2 ):0.001:pul( 2 );
prey = zeros( size( pylims ) );
sl = pll( [ 1 3 ] );
su = pul( [ 1 3 ] );
det_s = prod( su - sl );
spnodes = ( repmat( su, size( nodes, 1 ), 1 ) - repmat( sl, size( nodes, 1 ), 1 ) ) .* nodes + repmat( sl, size( nodes, 1 ), 1 );
parfor i = 1:length( prey )
    eval_sp = zeros( size( spnodes, 1 ), 3 );
    eval_sp( :, 1 ) = spnodes( :, 1 );
    eval_sp( :, 2 ) = pylims( i );
    eval_sp( :, 3 ) = spnodes( :, 2 );
    prey( i ) = sum( exp( dmsubs( lambdasc' * xhp, x, eval_sp' ) )' .* weights .* det_s ) * 1/prod( coord_scale );
end

pzlims = linspace( pll( 3 ), pul( 3 ), 1000 );
prez = zeros( size( pzlims ) );
sl = pll( [ 1 2 ] );
su = pul( [ 1 2 ] );
det_s = prod( su - sl );
spnodes = ( repmat( su, size( nodes, 1 ), 1 ) - repmat( sl, size( nodes, 1 ), 1 ) ) .* nodes + repmat( sl, size( nodes, 1 ), 1 );
parfor i = 1:length( prez )
    eval_sp = zeros( size( spnodes, 1 ), 3 );
    eval_sp( :, 1 ) = spnodes( :, 1 );
    eval_sp( :, 2 ) = spnodes( :, 2 );
    eval_sp( :, 3 ) = pzlims( i );
    prez( i ) = sum( exp( dmsubs( lambdasc' * xhp, x, eval_sp' ) )' .* weights .* det_s ) * 1/prod( coord_scale );
end

%% creating probability distributions in spherical space
[ TH, PHI, pthphi, R, pr ] = cart2sphPlot( lambdasc, coord_scale, x, xhp, 100 );

%% creating probability distributions on scaled data
xlims = ll( 1 ):0.1:ul( 1 );
px = zeros( size( xlims ) );
sl = ll( [ 2 3 ] );
su = ul( [ 2 3 ] );
det_s = prod( su - sl );
spnodes = ( repmat( su, size( nodes, 1 ), 1 ) - repmat( sl, size( nodes, 1 ), 1 ) ) .* nodes + repmat( sl, size( nodes, 1 ), 1 );
for i = 1:length( px )
    eval_sp = zeros( size( spnodes, 1 ), 3 );
    eval_sp( :, 1 ) = xlims( i );
    eval_sp( :, 2 ) = spnodes( :, 1 );
    eval_sp( :, 3 ) = spnodes( :, 2 );
    px( i ) = sum( exp( dmsubs( lambda' * xhp, x, eval_sp' ) )' .* weights .* det_s );
end

ylims = ll( 2 ):0.1:ul( 2 );
py = zeros( size( ylims ) );
sl = ll( [ 1 3 ]);
su = ul( [ 1 3 ]);
det_s = prod( su - sl );
spnodes = ( repmat( su, size( nodes, 1 ), 1 ) - repmat( sl, size( nodes, 1 ), 1 ) ) .* nodes + repmat( sl, size( nodes, 1 ), 1 );
for i = 1:length( py )
    eval_sp = zeros( size( spnodes, 1 ), 3 );
    eval_sp( :, 1 ) = spnodes( :, 1 );
    eval_sp( :, 2 ) = ylims( i );
    eval_sp( :, 3 ) = spnodes( :, 2 );
    py( i ) = sum( exp( dmsubs( lambda' * xhp, x, eval_sp' ) )' .* weights .* det_s );
end

zlims = ll( 3 ):0.01:ul( 3 );
pz = zeros( size( zlims ) );
sl = ll( [ 1 2 ]);
su = ul( [ 1 2 ]);
det_s = prod( su - sl );
spnodes = ( repmat( su, size( nodes, 1 ), 1 ) - repmat( sl, size( nodes, 1 ), 1 ) ) .* nodes + repmat( sl, size( nodes, 1 ), 1 );
parfor i = 1:length( pz )
    eval_sp = zeros( size( spnodes, 1 ), 3 );
    eval_sp( :, 1 ) = spnodes( :, 1 );
    eval_sp( :, 2 ) = spnodes( :, 2 );
    eval_sp( :, 3 ) = zlims( i );
    pz( i ) = sum( exp( dmsubs( lambda' * xhp, x, eval_sp' ) )' .* weights .* det_s );
end
%% plotting functionality 
[ histFreq, histXout ] = hist( predata( :, 1 ), 10 );
binWidth = histXout( 2 ) - histXout ( 1 );
figure; plot( histXout, histFreq/binWidth/sum( histFreq ), 'r' ); hold on;
plot( pxlims, prex, 'b' );

[ histFreq, histXout ] = hist( predata( :, 2 ), 10 );
binWidth = histXout( 2 ) - histXout( 1 );
figure; plot( histXout, histFreq/binWidth/sum( histFreq ), 'r' ); hold on;
plot( pylims, prey, 'b' );

[ histFreq, histXout ] = hist( predata( :, 3 ), 10 );
binWidth = histXout( 2 ) - histXout( 1 );
figure; plot( histXout, histFreq/binWidth/sum( histFreq ), 'r' ); hold on;
plot( pzlims, prez, 'b' );

[ X, Y, Z ] = sph2cart( TH, PHI, 1 );
figure; surf( X, Y, Z, pthphi );

[ histFreq, histXout ] = hist( sphdata( :, 3 ), 10 );
binWidth = histXout( 2 ) - histXout( 1 );
figure; plot( histXout, histFreq/binWidth/sum( histFreq ), 'r' ); hold on;
plot( R, pr, 'b' );

[ histFreq, histXout ] = hist( data( :, 3 ), 10 );
binWidth = histXout( 2 ) - histXout( 1 );
figure; plot( histXout, histFreq/binWidth/sum( histFreq ), 'r' ); hold on;
plot( zlims, pz, 'b' );
