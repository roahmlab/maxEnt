   function [ data, postdata, coord_scale ] = scaleDataPER( predata, per )
%
%  predata  -- struct with Talia's processed data
%  per      -- the percentile to choose, everyting is scaled to be between
%  [-1,1] 
%
% Create a single n-by-3 array of data points with the amount of scaling
% applied in each coordinate


% scaling_helper = range( predata );
coord_scale = prctile( abs( predata ), per );
% idx = all( abs( predata ) <= repmat( coord_scale, size( predata, 1 ), 1 ), 2 );
% predata = predata( idx, : );

postdata = predata;
coord_scale = coord_scale;
% coord_scale = scaling_helper/scale; 

% apply scaling 
for i = 1:size( predata, 2 )
    data( :, i ) = predata( :, i )/coord_scale( i );
end
