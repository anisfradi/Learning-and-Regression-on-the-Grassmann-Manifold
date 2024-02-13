
%
% INPUTS
% nIters = the number of iterations 
% obs = the set of time instances with corresponding object on the
% Grassmannian
% value_init = initialisation for position and velocity
% step = the step-size
%
% 
% OUTPUTS
% chain = a list of sampling for position and velocity
% acceptance = acceptance rate of MCMC 
% like = a list of log-likelihood values for each sample


function [chain, acceptance, like] = Grassmann_MCMC(nIters, obs, value_init,step)

chain = cell(nIters,1); like=[];
theta = value_init;
index = 0;
for iI = 1:nIters
    proposal = proposal(theta,step);
    propOld = Posterior(theta, obs);
    propNew = Posterior(proposal, obs);
    proba = exp(propNew - propOld);
    u=unifrnd(0, 1);
    if(u < proba)
        index = index + 1;
        theta=proposal;
    else
        index=index;
    end
    fprintf('Iteration %d\n', iI);
    chain{iI}=theta;
    like(iI)=Likelihood(chain{iI}, obs);
end

acceptance = index / nIters;

end


function prop = proposal(params,step)

pos = params.X1;
vel = params.X2;

% generate a sample on manifold
[n, p] = size(pos);
tan = normrnd(0, 0.001, [n, p]);
tan = (eye(n) - pos * pos') * tan;
[~, Ytmp, ~] = integrateForwardWithODE45(pos, tan, [0, 1]);
new_pos = Ytmp{end};
prop.X1 = new_pos;

[n, p] = size(vel);
new_vel = normrnd(vel, 0.001, [n, p]);
new_vel = (eye(n) - pos * pos') * new_vel;
[~, ~, ~, tan] = grgeo(pos, new_pos, step, 'v3', 'v2');
new_vel = ptEdelman(pos, tan, new_vel, step); 

prop.X1 = new_pos;
prop.X2 = new_vel;

end

function posterior = Posterior(params, obs)

posterior = Prior(params) + Likelihood(params, obs);

end

function prior = Prior(params)

pos = params.X1;
vel = params.X2;

%  posPrior = log(normpdf(pos, 0, 0.1));
%  velPrior = log(normpdf(vel, 0, 0.1));
% 
%  prior = sum(sum(posPrior)) + sum(sum(velPrior));


posPrior = -10*norm(pos, 'fro')^2/2;
velPrior = -10*norm(vel, 'fro')^2/2;

prior = posPrior + velPrior;

end

function likelihood = Likelihood(params, obs)

grParams = setRegressParams();
grParams.Ys = obs.shapes;
grParams.ts = obs.ts;
grParams.wi = ones(length(grParams.ts), 1);
[grParams.ts, idSort] = sort(grParams.ts);
grParams.Ys = grParams.Ys(idSort);
grParams.wi = grParams.wi(idSort);
grParams.ts_pre = grParams.ts;   
grParams.ts = grParams.ts - grParams.ts(1);
grParams.ts = grParams.ts ./ grParams.ts(end);
grParams.ts = min( max( round(grParams.ts / grParams.h) + 1, 1 ), grParams.nt+1 );

[~, X1s, X2s] = integrateForwardWithODE45(params.X1, params.X2, (0:grParams.nt)*grParams.h);
sigmaSq = estimateSigmaOfNoiseModel(X1s(grParams.ts), grParams.Ys);

energyV = grParams.alpha * trace( (X2s{1})' * X2s{1} );
energyS = 0;
for iI = 1:length(grParams.ts)
    [~, ~, ~, vTmp] = grgeo(X1s{grParams.ts(iI)}, grParams.Ys{iI}, 1, 'v3', 'v2');
    tmpDistSq = computeGeodesciDistanceWithSigma(sigmaSq{iI}, vTmp);
    energyS = energyS + tmpDistSq * grParams.wi(iI);
end
energyS = energyS / 2.0;  % here we compute energy / 2sigma^2
likelihood = - energyV - energyS;

end

