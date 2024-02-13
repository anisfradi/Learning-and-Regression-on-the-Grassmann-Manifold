% compute the distances between geodesic and measurements
% Input:
% Y0: the initial point on the manifold
% Y0dot: the initial velocity in the tangent space
% t0: the associated value of the initial point
% oma: the measurements
% t: the sequence of the associated values
% tStep: the time step for integrating forward
%
% Output:
% distSqrSum: the sum of square distances
% 
% Note: use ODE45 to integrateForward
% 
% Author: Yi Hong
% Date: 10/02/2013
%
function distSqrSum = sumDistOfPointsToGeodesicODE45(Y0, Y0dot, t0, oma, t, distSumMin, useODE45, weight )

if nargin <= 5
    distSumMin = -1;
    useODE45 = true;
    weight = ones(size(t));
end

maxT = max(t);
minT = min(t);
hStep = (maxT - minT) / 100.0;

% integrate to the boundary of time range
if t0 == minT
    if useODE45
        %tspan = linspace( t0, max(t), 50 );
        tspan = t0:hStep:maxT;
        if tspan(end) < maxT
            tspan(end+1) = maxT;
        end
        %tspan
        [tForward, YForward] = integrateForwardWithODE45(Y0, Y0dot, tspan);
    else
        [YForward, ~, tForward] = integrateForwardToGivenTime(Y0, Y0dot, max(t)-t0, (max(t)-t0)/50.0 );
        tForward = tForward + t0;
    end
    tBackward(1) = t0;
    YBackward{1} = Y0;
elseif t0 == maxT
    if useODE45
        %tspan = linspace( t0, min(t), 50 );
        tspan = t0:-hStep:minT;
        if tspan(end) > minT
            tspan(end+1) = minT;
        end
        [tBackward, YBackward] = integrateForwardWithODE45(Y0, Y0dot, tspan);
    else
        [YBackward, ~, tBackward] = integrateForwardToGivenTime(Y0, Y0dot, min(t)-t0, (t0-min(t))/50.0);
        tBackward = tBackward + t0;
    end
    tForward(1) = t0;
    YForward{1} = Y0;
else
    if useODE45
        %tspan = linspace( t0, max(t), 50 );
        tspan = t0:hStep:maxT;
        if tspan(end) < maxT
            tspan(end+1) = maxT;
        end
        [tForward, YForward] = integrateForwardWithODE45(Y0, Y0dot, tspan);
        %tspan = linspace( t0, min(t), 50 );
        tspan = t0:-hStep:minT;
        if tspan(end) > minT
            tspan(end+1) = minT;
        end
        [tBackward, YBackward] = integrateForwardWithODE45(Y0, Y0dot, tspan);
    else
        [YForward, ~, tForward] = integrateForwardToGivenTime(Y0, Y0dot, max(t)-t0, (max(t)-t0)/50.0 );
        tForward = tForward + t0;
        [YBackward, ~, tBackward] = integrateForwardToGivenTime(Y0, Y0dot, min(t)-t0, (t0-min(t))/50.0);
        tBackward = tBackward + t0;
    end
end
        
% compute the sum of the square distances
distSqrSum = 0;
for id = 1:length(oma)
    if t(id) >= t0
        % find t(id) index in the forward results
        [~, tIndex] = min( abs(tForward - t(id)) );
        oma_t = YForward{tIndex};
    else
        % find t(id) index in the backward results
        [~, tIndex] = min( abs(tBackward - t(id)) );
        oma_t = YBackward{tIndex};
    end
    distSqrSum = distSqrSum + grarc( oma{id}, oma_t ) * weight(id);
    % terminate the computation if the sum is over the current one 
    if distSumMin > 0 && distSumMin < distSqrSum
        break;
    end
end
