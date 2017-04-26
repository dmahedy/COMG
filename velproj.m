function [ p, n ] = velproj( u, v, uvb )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     [ uvb, uvstats ] = robustfit(u,v);

    p = NaN( size( u, 1 ), size( u, 2 ) );
    n = NaN( size( u, 1 ), size( u, 2 ) );
    for imx = 1:size(u,2)
        domvc = [ (max(u(:,imx)) - 0), ...
            ((uvb(2) * max(u(:,imx)) + uvb(1) ) - uvb(1) ) ];
        domvc = domvc / ( sqrt( domvc( 1 )^2 + domvc( 2 )^2 ) );
        pvect = ( dot( repmat( domvc, length( u(:,imx) ), 1 ), ...
            [ u( :, imx ), v( :, imx ) ], 2 ) ) .* ...
            repmat( domvc, length( u( :, imx ) ), 1 );
        nvect = [ u( :, imx ), v( :, imx ) ] - pvect;
        pt = sqrt( pvect( :, 1 ).^2 + pvect( :, 2 ).^2 );
        nt = sqrt( nvect( :, 1 ).^2 + nvect( :, 2 ).^2 );
        for npsg = 1:length( pvect )
            if all( pvect( npsg, 1 ) < 0 && pvect( npsg, 2 ) < 0 ) || ...
                all( pvect( npsg, 1 ) < 0 && abs( pvect( npsg, 1 ) ) > abs( pvect( npsg, 2 ) ) ) || ...
                all( pvect( npsg, 2 ) < 0 && abs( pvect( npsg, 2 ) ) > abs( pvect( npsg, 1 ) ) )
                    pt( npsg, 1 ) = -pt( npsg, 1 );
            end
            if all( nvect( npsg, 1 ) < 0 && nvect( npsg, 2 ) < 0 ) || ...
                all( nvect( npsg, 1 ) < 0 && abs( nvect( npsg, 1 ) ) > abs( nvect( npsg, 2 ) ) ) || ...
                all( nvect( npsg, 2 ) < 0 && abs( nvect( npsg, 2 ) ) > abs( nvect( npsg, 1 ) ) )
                    nt( npsg, 1 ) = -nt( npsg, 1 );
            end
        end
        
        p( :, imx ) = pt;
        n( :, imx ) = nt;
    end

end

