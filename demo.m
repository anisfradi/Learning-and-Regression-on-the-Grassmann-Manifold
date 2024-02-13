
clear all;
close all;
clc

%% set options of Gr(n,p)
N = 21; % sample size
n = 3; p = 1; % matrix n x p
G = grassmannfactory(n,p); % geometric tools of the Grassmannian Gr(n,p)
time = linspace(0,1,N); % time instances between 0 and 1


%% generate N random points on Gr(n,p) following the geodesic path
norm_vel=pi/2; % injectivity radius

while norm_vel >= pi/2
p0_true=G.rand(); % initial position
sim=normrnd(0, 1, [n,p]); 
v0_true=(eye(n) - p0_true * p0_true') * sim; % initial velocity
norm_vel=norm(v0_true,'fro');
end

p1_true=G.exp(p0_true,v0_true); % end position

var_noise = 0.1; % noise variance

for i=1:N
    Q_true{i}=G.exp(p0_true,time(i)*v0_true); % clean data
end

Q=addNoiseToData(Q_true, var_noise); % add noise to clean data


%% plot simulated data on Gr(n=3,p=1)
plotSimulationsGr31

                                        
%% run MCMC on the Grassmannian 
obs.shapes = Q;
obs.ts = time;

% Initialization
setup_initialization
[p_init, v_init] = initializeRegression( params.Ys, params.ts, params.alpha, params.sigmaSqr, params.useODE45, params.useRansac, params.wi );
value_init.X1=p_init; value_init.X2=v_init; 

iter=1e3; % number of iterations for Monte Carlo sampling
[chain, acceptance, like] = Grassmann_MCMC(iter, obs, value_init,1);


%% Markov chain of position and velocity
chain_pos=cell(iter,1); chain_poss=[]; chain_vel=[];
for i=1:iter
chain_pos{i}=chain{i}.X1;
chain_poss(:,i)=chain{i}.X1;
chain_vel(:,i)=chain{i}.X2;
end

%% plot one example of Markov Chain position (top) and velocity and their density estimations (bottom)
exp_pos=reshape(chain_poss(1,:),[iter,1]);
exp_vel=reshape(chain_vel(1,:),[iter,1]);
[pdf_pos,support_pos] = ksdensity(exp_pos); 
[pdf_vel,support_vel] = ksdensity(exp_vel); 
plotMarkovChain_Gr31

%% find mean values of Markov Chain position and velocity to be the initial condition estimations
mean_pos=grfmean(chain_pos, 1e-5); % Fréchet mean
mean_vel=mean(chain_vel,2); mean_vel=reshape(mean_vel,[n,p]); % Euclidean mean

%% Difference in initial conditions
diff_pos = sqrt(grarc(mean_pos, p0_true));
sprintf(['difference norm in position is: ', num2str(diff_pos)])
diff_vel = norm(mean_vel-v0_true, 'fro');
sprintf(['difference norm in velocity is: ', num2str(diff_vel)])

%% plot the cost function: the minus of log-likelihood
figure('DefaultAxesFontSize',25,'DefaultAxesFontWeight','bold');
plot(1:iter,-like,'-b','Linewidth',1)
xlabel('iterations')
set(gcf,'color','w'); % Set background color to white
set(gca,'box','off'); % Turn off the axes box

%% plot fitting on Gr(n=3,p=1)
plotFittingGr31


%% plot 2D Grassmannian distance
figure('DefaultAxesFontSize',25,'DefaultAxesFontWeight','bold');
plot2DGrass(fitQ, Q, time)
xlabel('time')
set(gcf,'color','w'); % Set background color to white
set(gca,'box','off'); % Turn off the axes box

%% compute the R-squared coefficient
R2 = 1 - computeVar(Q,fitQ, 3) / computeVar(Q, Q, 1);
sprintf(['R-squared coefficient is: ', num2str(R2)])

%% compute the mean squared error (MSE)
MSE = mse(G,Q,mean_pos,mean_vel,time);
sprintf(['mean squared error is: ', num2str(MSE)])

%% compute the data-to-noise-ratio (DNR)
SNR = computeVar(Q, Q, 1) / computeVar(Q,fitQ, 3);
sprintf(['DNR is: ', num2str(SNR)])

