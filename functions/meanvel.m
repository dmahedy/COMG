function [ vnint, vpint, tngdif, tpsdif, vnmean, vpmean, ...
    npkpts, ppkpts, vnpeak, vppeak ] = meanvel( ts, velmat )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    vnint = NaN( 1, size( velmat, 2 ) );
    vpint = NaN( 1, size( velmat, 2 ) );
    
%     tngdif = NaN( 1, size( velmat, 2 ) );
%     tpsdif = NaN( 1, size( velmat, 2 ) );
    
    vnmean = NaN( 1, size( velmat, 2 ) );
    vpmean = NaN( 1, size( velmat, 2 ) );
    
    npkpts = cell( 1, size( velmat, 2 ) );
    ppkpts = cell( 1, size( velmat, 2 ) );
    
    vnpeak = NaN( 1, size( velmat, 2 ) );
    vppeak = NaN( 1, size( velmat, 2 ) );
    
    for lyr = 1:size( velmat, 2 )
        ttmp = ts;
        veltmp = velmat( :, lyr );
        
        figure;
        plot( ttmp, veltmp );
        hold on;
        plot( [ ttmp( find( isnan( veltmp ) )-1 ); ...
            ttmp( find( isnan( veltmp ) )+1 ) ], ...
            [ veltmp( find( isnan( veltmp ) )-1 ); ...
            veltmp( find( isnan( veltmp ) )+1 ) ], 'ob' );
        
        ttmp2 = ttmp; veltmp2 = veltmp;
        ttmp2( isnan( veltmp2 ) ) = [];
        veltmp2( isnan( veltmp2 ) ) = [];

        [ t0, v0, i, ~ ] = intersections( ttmp2, veltmp2, ...
            [min(ttmp2) max(ttmp2)], [0 0], 'true' );
        
        t0( i( mod(i,1) == 0 ) ) = [];
        v0( i( mod(i,1) == 0 ) ) = [];
        i( i( mod(i,1) == 0 ) ) = [];
        
        t0( isnan( i ) ) = [];
        v0( isnan( i ) ) = [];
        i( isnan( i ) ) = [];
        
        plot( t0, v0, 'sk' );
        
        if sign( veltmp( floor( i(1) ) ) ) == 1 && ...
                mod( numel( t0 ), 2 ) == 0
            tngdif = t0( 2:2:end ) - t0( 1:2:end-1 );
            tpsdif = t0( 3:2:end-1 ) - t0( 2:2:end-2 );
        elseif sign( veltmp( floor( i(1) ) ) ) == 1 && ...
                mod( numel( t0 ), 2 ) == 1
            tngdif = t0( 2:2:end-1 ) - t0( 1:2:end-2 );
            tpsdif = t0( 3:2:end ) - t0( 2:2:end-1 );
        elseif sign( veltmp( floor( i(1) ) ) ) == -1 && ...
                mod( numel( t0 ), 2 ) == 0
            tngdif = t0( 3:2:end-1 ) - t0( 2:2:end-2 );
            tpsdif = t0( 2:2:end ) - t0( 1:2:end-1 );
        elseif sign( veltmp( floor( i(1) ) ) ) == -1 && ...
                mod( numel( t0 ), 2 ) == 1
            tngdif = t0( 3:2:end ) - t0( 2:2:end-1 );
            tpsdif = t0( 2:2:end-1 ) - t0( 1:2:end-2 );
        end
         
        [ tspos, idtp ] = sort( [ ttmp( sign( veltmp ) == 1 ); t0 ] );
        velpos = [ veltmp( sign( veltmp ) == 1 ); v0 ]; 
        velpos = velpos( idtp );
        
        [ tsneg, idtn ] = sort( [ ttmp( sign( veltmp ) == -1 ); t0 ] );
        velneg = [ veltmp( sign( veltmp ) == -1 ); v0 ]; 
        velneg = velneg( idtn );
        
        plot( tspos, velpos, '-g', 'LineWidth', 1 );
        plot( tsneg, velneg, '-r', 'LineWidth', 1 );
        
        vnint( 1, lyr ) = trapz( tsneg, velneg );
        vpint( 1, lyr ) = trapz( tspos, velpos );
        
        vnmean( 1, lyr ) = vnint( 1, lyr ) / sum( tngdif );
        vpmean( 1, lyr ) = vpint( 1, lyr ) / sum( tpsdif );
        
        %% Peak Currents
        
        tavdif = median( [ tngdif; tpsdif ] );
        
        [vpkps,tspkps] = findpeaks( veltmp2( ceil( i(1) ):floor( i(end) ) ), ...
            ttmp2( ceil( i(1) ):floor( i(end) ) ), ...
            'MinPeakDistance', tavdif );
        [vpkng,tspkng] = findpeaks( -veltmp2( ceil( i(1) ):floor( i(end) ) ), ...
            ttmp2( ceil( i(1) ):floor( i(end) ) ), ...
            'MinPeakDistance', tavdif );
        vpkng = vpkng * -1;
        tspkps( vpkps < 0 ) = []; vpkps( vpkps < 0 ) = [];
        tspkng( vpkng > 0 ) = []; vpkng( vpkng > 0 ) = [];
        
        plot( tspkps, vpkps, '*r', 'LineWidth', 2 );
        plot( tspkng, vpkng, 'og', 'LineWidth', 2 );
        
        npkpts( 1, lyr ) = { [ tspkps, vpkps ] };
        ppkpts( 1, lyr ) = { [ tspkng, vpkng ] };
        
        vnpeak( 1, lyr ) = mean( vpkng );
        vppeak( 1, lyr ) = mean( vpkps );
        
    end

end

