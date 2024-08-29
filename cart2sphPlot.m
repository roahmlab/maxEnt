function [ TH, PHI, pthphi, R, pr ] = cart2sphPlot( lambdasc, coord_scale, x, xhp, RES )

R = linspace( 0, sqrt( coord_scale * coord_scale' ), RES );
pr = zeros( size( R ) );
parfor i = 1:length( R )
    helperphi1 = asin( min( coord_scale( 3 )/R( i ), 1 ) );
    helperphi2 = pi - asin( min( coord_scale( 3 )/R( i ), 1 ) );
    helperphi3 = pi - asin( max( -coord_scale( 3 )/R( i ), -1 ) );
    helperphi4 = 2 * pi + asin( max( -coord_scale( 3 )/R( i ), -1 ) );
    
    
    PHI = [ linspace( 0, helperphi1, ceil( RES/3 ) )  linspace( helperphi2, helperphi3, ceil( RES/3 ) ) ...
        linspace( helperphi4, 2 * pi, ceil( RES/3 ) ) ];
    intphi = zeros( size( PHI ) );
    for j = 1:length( PHI )
        helperx1 = acos( min( coord_scale( 1 )/( R( i ) * cos( PHI( j ) ) ), 1 ) );
        helperx2 = acos( max( -coord_scale( 1 )/( R( i ) * cos( PHI( j ) ) ), -1 ) );
        helperx3 = pi + acos( min( coord_scale( 1 )/( R( i ) * cos( PHI( j ) ) ), 1 ) );
        helperx4 = pi + acos( max( -coord_scale( 1 )/( R( i ) * cos( PHI( j ) ) ), -1 ) );
        helpery1 = asin( min( coord_scale( 2 )/( R( i ) * cos( PHI( j ) ) ), 1 ) );
        helpery2 = pi - asin( min( coord_scale( 2 )/( R( i ) * cos( PHI( j ) ) ), 1 ) );
        helpery3 = pi - asin( max( -coord_scale( 2 )/( R( i ) * cos( PHI( j ) ) ), -1 ) );
        helpery4 = 2 * pi + asin( max( -coord_scale( 2 )/( R( i ) * cos( PHI( j ) ) ), -1 ) );
        
        if( helperx1 < helpery1 )
            TH = linspace( helperx1, helpery1, RES/4 );
            X = R( i ) * cos( PHI( j ) ) * cos( TH );
            Y = R( i ) * cos( PHI( j ) ) * sin( TH );
            Z = repmat( R( i ) * sin( PHI( j ) ), size( TH ) );
            helper = exp( dmsubs( lambdasc' * xhp, x, [ X; Y; Z ] ) ) * 1/prod( coord_scale ) .* ...
                ( ( X ).^2 + ( Y ).^2 + ( Z ).^2 ) .* ...
                cos( atan2( ( Z ), sqrt( ( X ).^2 + ( Y ).^2 ) ) );
%                 ( ( X * coord_scale( 1 ) ).^2 + ( Y* coord_scale( 2 ) ).^2 + ( Z * coord_scale( 3 ) ).^2 ) .* ...
%                 cos( atan2( ( Z * coord_scale( 3 ) ), sqrt( ( X * coord_scale( 1 ) ).^2 + ( Y* coord_scale( 2 ) ).^2 ) ) );
            intphi( j ) = intphi( j ) + sum( helper ) * ( TH( 2 ) - TH( 1 ) );
        end
        if( helpery2 < helperx2 )
            TH = linspace( helpery2, helperx2, RES/4 );
            X = R( i ) * cos( PHI( j ) ) * cos( TH );
            Y = R( i ) * cos( PHI( j ) ) * sin( TH );
            Z = repmat( R( i ) * sin( PHI( j ) ), size( TH ) );
            helper = exp( dmsubs( lambdasc' * xhp, x, [ X; Y; Z ] ) ) * 1/prod( coord_scale ) .* ...
                ( ( X ).^2 + ( Y ).^2 + ( Z ).^2 ) .* ...
                cos( atan2( ( Z ), sqrt( ( X ).^2 + ( Y ).^2 ) ) );
%                 ( ( X * coord_scale( 1 ) ).^2 + ( Y* coord_scale( 2 ) ).^2 + ( Z * coord_scale( 3 ) ).^2 ) .* ...
%                 cos( atan2( ( Z * coord_scale( 3 ) ), sqrt( ( X * coord_scale( 1 ) ).^2 + ( Y* coord_scale( 2 ) ).^2 ) ) );
            intphi( j ) = intphi( j ) + sum( helper ) * ( TH( 2 ) - TH( 1 ) );
        end
        if ( helperx3 < helpery3 )
            TH = linspace( helperx3, helpery3, RES/4 );
            X = R( i ) * cos( PHI( j ) ) * cos( TH );
            Y = R( i ) * cos( PHI( j ) ) * sin( TH );
            Z = repmat( R( i ) * sin( PHI( j ) ), size( TH ) );
            helper = exp( dmsubs( lambdasc' * xhp, x, [ X; Y; Z ] ) ) * 1/prod( coord_scale ) .* ...
                ( ( X ).^2 + ( Y ).^2 + ( Z ).^2 ) .* ...
                cos( atan2( ( Z ), sqrt( ( X ).^2 + ( Y ).^2 ) ) );
%                 ( ( X * coord_scale( 1 ) ).^2 + ( Y* coord_scale( 2 ) ).^2 + ( Z * coord_scale( 3 ) ).^2 ) .* ...
%                 cos( atan2( ( Z * coord_scale( 3 ) ), sqrt( ( X * coord_scale( 1 ) ).^2 + ( Y* coord_scale( 2 ) ).^2 ) ) );
            intphi( j ) = intphi( j ) + sum( helper ) * ( TH( 2 ) - TH( 1 ) );
        end
        if( helpery4 < helperx4 )
            TH = linspace( helpery4, helperx4, RES/4 );
            X = R( i ) * cos( PHI( j ) ) * cos( TH );
            Y = R( i ) * cos( PHI( j ) ) * sin( TH );
            Z = repmat( R( i ) * sin( PHI( j ) ), size( TH ) );
            helper = exp( dmsubs( lambdasc' * xhp, x, [ X; Y; Z ] ) ) * 1/prod( coord_scale ) .* ...
                ( ( X ).^2 + ( Y ).^2 + ( Z ).^2 ) .* ...
                cos( atan2( ( Z ), sqrt( ( X ).^2 + ( Y ).^2 ) ) );
%                 ( ( X * coord_scale( 1 ) ).^2 + ( Y* coord_scale( 2 ) ).^2 + ( Z * coord_scale( 3 ) ).^2 ) .* ...
%                 cos( atan2( ( Z * coord_scale( 3 ) ), sqrt( ( X * coord_scale( 1 ) ).^2 + ( Y* coord_scale( 2 ) ).^2 ) ) );
            intphi( j ) = intphi( j ) + sum( helper ) * ( TH( 2 ) - TH( 1 ) );
        end
    end
    
    pr( i ) = sum( intphi( 1:ceil( RES/3 ) ) ) * ( PHI( 2 ) - PHI( 1 ) ) +  ...
        sum( intphi( ( ceil( RES/3 ) + 1 ):( 2 * ceil( RES/3 ) ) ) ) * ( PHI( ceil( RES/3 ) + 2 ) - PHI( ceil( RES/3 ) + 1 ) ) + ...
        sum( intphi( ( 2 * ceil( RES/3 ) + 1 ):end ) ) * ( PHI( ( 2 * ceil( RES/3 ) + 2 ) ) - PHI( ( 2 * ceil( RES/3 ) + 1 ) ) );
end

thlims = linspace( 0, 2 * pi + 0.1, RES );
philims = linspace( -pi/2, pi/2, RES );
pthphi = zeros( length( thlims ), length( philims ) );
TH = zeros( length( thlims ), length( philims ) );
PHI = zeros( length( thlims ), length( philims ) );
for i = 1:length( PHI )
    for j = 1:length( TH )
        R3 = abs( coord_scale( 3 )/sin( philims( i ) ) );
        R2 = abs( coord_scale( 2 )/( cos( philims( i ) ) * sin( thlims( j ) ) ) );
        R1 = abs( coord_scale( 1 )/( cos( philims( i ) ) * cos( thlims( j ) ) ) );
        R = linspace( 0, min( [ R3 R2 R1 ] ) );
        
        X = R * cos( philims( i ) ) * cos( thlims( j ) );
        Y = R * cos( philims( i ) ) * sin( thlims( j ) );
        Z = R * sin( philims( i ) );
        helper = exp( dmsubs( lambdasc' * xhp, x, [ X; Y; Z ] ) ) * 1/prod( coord_scale ) .* ...
            ( ( X * coord_scale( 1 ) ).^2 + ( Y* coord_scale( 2 ) ).^2 + ( Z * coord_scale( 3 ) ).^2 ) .* ...
            cos( atan2( ( Z * coord_scale( 3 ) ), sqrt( ( X * coord_scale( 1 ) ).^2 + ( Y* coord_scale( 2 ) ).^2 ) ) );
        pthphi( j, i ) = sum( helper ) * ( R( 2 ) - R( 1 ) );
        TH( j, i ) = thlims( j );
        PHI( j, i ) = philims( i );
    end
end

R = linspace( 0, sqrt( coord_scale * coord_scale' ), RES );
