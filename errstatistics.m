function [ Rsq, RMSE, SI, bias, bspct, resid ] = errstatistics( t1, p1, t2, p2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    t1( isnan( p1 ) ) = []; t2( isnan( p2 ) ) = [];
    p1( isnan( p1 ) ) = []; p2( isnan( p2 ) ) = [];

    [b, bstat] = robustfit( p1( ismember( t1, t2 ) ), ...
        p2( ismember( t2, t1 ) ) );
    Rsq = corr( p2( ismember( t2, t1 ) ), ...
        b(1) + b(2) * p1( ismember( t1, t2 ) ))^2;
    RMSE = bstat.s;
    SI = bstat.s / ...
        mean( p1( ismember( t1, t2 ) ) );
    bias = mean( bstat.resid );
    bspct = 100*( sum(p2( ismember( t2, t1 ) )) - ...
        sum(p1( ismember( t1, t2 ) )) ) ./ ...
        sum(p1( ismember( t1, t2 ) ));
    resid = bstat.resid;

end

