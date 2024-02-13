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

function [Y0, Y0dot, t0, distMin] = grpairSearchODE45( oma, t, nFlag, alpha, sigmaSqr, useODE45, useRansac, weight )

%strTmp = sprintf('Compute geodesic of %d objects using pair-search (geodesic defined by two points)', length(t));
%disp(strTmp);

nRandPoints = 2;        % the number of points sampled for regression
nPoints = length(t);    % the total number of points

if (~useRansac || useRansac && (nPoints < 40)) 
    % directly compute the geodesic using pair-wise searching
    [Y0, Y0dot, t0, distMin] = subPairSearch( oma, t, nFlag, oma, t, alpha, sigmaSqr, useODE45, weight );
else
    %disp( 'Using random sample for pair-wise searching...' );
    iterMax = 100; 
    iterMaxLoop = 10;
    distMin = -1;
    thresholdTmp = (max(t) - min(t))*1e-5;
    for iter = 1:iterMax
        rng('shuffle');
        idPerm = randperm(nPoints);
        idLoop = 1;
        while idLoop <= iterMaxLoop && abs( t(idPerm(1)) - t(idPerm(2)) ) < thresholdTmp
            idPerm = randperm(nPoints);
            idLoop = idLoop + 1;
        end
        %strTmp = sprintf( 'Iter: %d, t0: %f, t1: %f', iter, t(idPerm(1)), t(idPerm(2)) );
        %disp(strTmp);
        [Y0tmp, Y0dottmp, t0tmp, disttmp] = subPairSearch( oma(idPerm(1:nRandPoints)), ...
            t(idPerm(1:nRandPoints)), nFlag, oma, t, alpha, sigmaSqr, useODE45, weight );
        % the selected pair is too close, reselect
        %if disttmp == -1
        %    iter = iter - 1;
        %    continue;
        %end
        if distMin == -1 || ( disttmp > 0 && distMin > disttmp )
            Y0 = Y0tmp;
            Y0dot = Y0dottmp;
            t0 = t0tmp;
            distMin = disttmp;
        end        
    end
end

% integrate backward to get the initial conditions at min(t)
minT = min(t);
if t0 ~= minT
    if useODE45
        [~, Ytmp, Ydottmp] = integrateForwardWithODE45( Y0, -Y0dot, [0 (t0-minT)/2.0 t0-minT] );
    else
        [Ytmp, Ydottmp] = integrateForwardToGivenTime( Y0, -Y0dot, t0-minT, (t0-minT)/50.0 );
    end
    Y0 = Ytmp{end};
    Y0dot = -Ydottmp{end};
    t0 = minT;
end

    
end

% direct pair-wise searching
function [Y0, Y0dot, t0, distMin] = subPairSearch( oma, t, nFlag, omaAll, tAll, alpha, sigmaSqr, useODE45, weight )

    distMin = -1;
    for iI = 1:length(t)-1
        for iJ = length(t):-1:iI+1
            
            % do not use two points that are too close to each other
            if abs( t(iJ) - t(iI) ) < 1e-3
                continue;
            end
            [~, ~, ~, velInit] = grgeo( oma{iI}, oma{iJ}, 1, 'v3', 'v2' );
            % adjust the velocity to make it shooting from ti to tj
            velInit = velInit / (t(iJ) - t(iI));

            % compute the sum of the square distances
            if nFlag == 1   % compute the distance between measurements and the corresponding points on the geodesic
                distSum = sumDistOfPointsToGeodesicODE45(oma{iI}, velInit, t(iI), omaAll, tAll, -1, useODE45, weight);
            elseif nFlag == 2     % compute the distance of measurements projected to the geodesic
                distSum = sumDistOfPointsProjectToGeodesicODE45(oma{iI}, velInit, t(iI), omaAll, tAll, -1);
            else
                error('Not supported distance calculation');
            end
            distSum = alpha * trace(velInit'*velInit) + distSum ./ sigmaSqr;

            if distMin == -1 || distMin > distSum
                distMin = distSum;
                Y0 = oma{iI};
                Y0dot = velInit;
                t0 = t(iI);
                %t1 = t(iJ);
            end
        end
    end
end

