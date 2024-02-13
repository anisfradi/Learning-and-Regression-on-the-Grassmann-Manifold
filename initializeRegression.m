% geodesic defined by two points
% input:
% oma: the sequence of obserbility matrices
% t:   the associated value for each LDS
% tStep: the time step for integrating forward
% nFlag: 1 -- distance at given time points, 2 -- projection distance
% output:
% Y0: the initial point on the manifold
% Y0dot: the initial velocity for shooting
% t0: the associated value for the initial point
%
% Note: the distance of method 2 is measured by projection distances
% 
% Author: Yi Hong
% Date: 10/02/2013

function [Y0, Y0dot] = initializeRegression( oma, t, alpha, sigmaSqr, useODE45, useRansac, weight )

[Y0, Y0dot, t0] = grpairSearchODE45( oma, t, 1, alpha, sigmaSqr, useODE45, useRansac, weight );
% integrate backward to get the initial conditions at min(t)
minT = min(t);
maxT = max(t);
if t0 > minT
    if useODE45
        [~, Ytmp, Ydottmp] = integrateForwardWithODE45( Y0, -Y0dot, [0 (t0-minT)/2.0 t0-minT] );
    else
        [Ytmp, Ydottmp] = integrateForwardToGivenTime( Y0, -Y0dot, t0-minT, (t0-minT)/50.0 );
    end
    Y0 = Ytmp{end};
    Y0dot = -Ydottmp{end};
elseif t0 < minT
    if params.useODE45
        [~, Ytmp, Ydottmp] = integrateForwardWithODE45( Y0, Y0dot, [0 (minT-t0)/2.0 minT-t0] );
    else
        [Ytmp, Ydottmp] = integrateForwardToGivenTime( Y0, Y0dot, minT-t0, (minT-t0)/50.0 );
    end
    Y0 = Ytmp{end};
    Y0dot = Ydottmp{end};
end
Y0dot = Y0dot * (maxT-minT);

