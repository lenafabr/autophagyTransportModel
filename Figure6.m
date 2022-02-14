%% script to reproduce plots in figure 6
% requires running Figure5.m to generate a savefile with fitted parameters
% requires running agent based simulations with fitted parameters using
% code found at https://github.com/lenafabr/particleDynamics1D
% required files: fitparams_onedone_nobranch.mat, 
clear;clc;
%% add scripts to path
addpath('./tools/')
%% set up options for calculation
options = struct();

options.nx = 1000; % number of points on the calculation grid

options.solvePfuse = true;
options.dotipfusion = true; % fuse particles in the distal tip
options.onedone = true; % allow AVs to only fuse once

numBranch = 0; % number of branches
xvals = linspace(0,1,options.nx);
dx = xvals(2)-xvals(1);
%% load fitted parameters
load('./data/fitparams_onedone_nobranch.mat','fitparams','params','dens')
%% run a sweep over kpy and pf, fixing kye to the optimal value
% calculate the fraction of AVs fused for each parameter combination
params.kye = fitparams(3);

pflist = linspace(0.001,0.015,41);
kpylist = logspace(log10(1),log10(200),42);

Ffdist = zeros(length(pflist),length(kpylist));
Ydist = Ffdist;
YPrat = Ffdist;
allparams = cell(length(pflist),length(kpylist));
alldensities = allparams;

for pc = 1:length(pflist)
  for kc = 1:length(kpylist)
    
    params.kpy = kpylist(kc); % lysosome leaving rate
    params.pf = pflist(pc);
    params.pfra = params.pf; params.pfrr = params.pf; params.pfrs = params.pf;
    params.pfaa = params.pf; params.pfar = params.pf; params.pfas = params.pf;
    params.pfsa = params.pf; params.pfsr = params.pf;
    params.pft = params.pf; % fusion probability within tip
   
    allparams{pc,kc} = params;
    
    [dens,success] = onedone_branches_solver(params,options);
    
    alldensities{pc,kc} = dens;
    xvals = dens.xvals;
    
    % AV distributions
    Btot = dens.Ba+dens.Bs+dens.Br;
    Ptot = dens.R +Btot;
    Stot = Btot;
    
    Ytot = dens.Ya+dens.Yr+dens.Ysa+dens.Ysr+dens.Yuw+dens.Yus;
    
    Butot = dens.Bau+dens.Bsu+dens.Bru;
    Putot = dens.Ru+Butot;
    Ftot = Ptot-Putot;
    YFtot = Ytot+Ftot; % LAMP1 puncta;
    
    % distal fraction fused
    distind = find(xvals<xd);
    
    Ptotdist = (sum(Ptot(distind))-0.5*(Ptot(distind(1))+Ptot(distind(end))))*dx;
    Fdist =  (sum(Ftot(distind))-0.5*(Ftot(distind(1))+Ftot(distind(end))))*dx;
    Ffdist(pc,kc) = (Fdist)./(Ptotdist);
    
    % distal LC3 density
    Pdist = (Ptotdist)/xd/params.Lreal;
    % distal LAMP1 density
    Ytotdist = (sum(Ytot(distind))-0.5*(Ytot(distind(1))+Ytot(distind(end))))*dx;
    YFtotdist = (sum(YFtot(distind))-0.5*(YFtot(distind(1))+YFtot(distind(end))))*dx;
    Ydist(pc,kc) = (YFtotdist+dens.Yt0)/xd/params.Lreal;
    
    % proximal ratio of lysosome:AV density
    proxind = find(xvals>(1-xd));
    YPrat(pc,kc) = sum(YFtot(proxind))/sum(Ptot(proxind));
    
    disp([pc kc params.kpy params.pfra Ffdist(pc,kc) Ydist(pc,kc) YPrat(pc,kc)])
  end
  save('./data/onedone_nobranch_scanpfkpy_kyefit.mat')
end
%% interpolate to get fraction fused as a function of lysosome density and for each pf

rhoylist = linspace(0.03,0.12,20);
kpyfitrhoy = zeros(length(pflist),length(rhoylist));
Fffitrhoy = kpyfitrhoy;
pffitFf = zeros(1,length(rhoylist));

Ffwant = 0.49; % desired fraction fused.
for rc = 1:length(rhoylist)
  rhoy = rhoylist(rc);
  try
    for pc = 1:length(pflist)
      % at this pf, interpolate for the kpy that gives the desired rhoy
      kpyfitrhoy(pc,rc) = interp1(Ydist(pc,:),log10(kpylist),rhoy);
      Fffitrhoy(pc,rc) = interp1(Ydist(pc,:),Ffdist(pc,:),rhoy);
    end
    
    % get the pf that gives the desired fraction fused
    ind = find(~isnan(Fffitrhoy(:,rc)));
    pffitFf(rc) = interp1(Fffitrhoy(ind,rc),log10(pflist(ind)),Ffwant);
    pffitFf(rc) = 10.^pffitFf(rc);
  end
end
%% plot fraction fused vs rhoy and pf
% figure 6 panel a
figure('Name','Figure 6(a)','NumberTitle','off')
pcolor(pflist,rhoylist,Fffitrhoy')
shading flat
set(gca,'XScale','linear','YScale','linear')
hold all
semilogx(pffitFf,rhoylist,'k','LineWidth',2)
plot(pflist,0.068*ones(size(pflist)),'r--')
hold off
cb = colorbar;
xlabel('fusion probability p_f')
ylabel('distal lysosome density \rho_y (\mum^{-1})')
title(sprintf('Lysosome exit rate k_e^y = %0.2f min^{-1}',params.kye*params.vpreal/params.Lreal*60))
xlim([3,14]*1e-3)
ax = gca;
ax.XTick = (3:3:12)*1e-3;
ax.YTick = (0.04:0.02:0.12);
caxis([0,0.8])
cb.Ticks = (0:0.2:0.8);
ylabel(cb,'Fraction of fused AVs f_f')
plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28)
%% for the fitted pf, for a range of kpy, fit the kye to get the
% right distal density; calculate proximal ratio
rhoywant = 0.068;

params.pf = pflist(pc);
params.pfra = params.pf; params.pfrr = params.pf; params.pfrs = params.pf;
params.pfaa = params.pf; params.pfar = params.pf; params.pfas = params.pf;
params.pfsa = params.pf; params.pfsr = params.pf;
params.pft = params.pf; % fusion probability within tip

kpylist2 = linspace(1,25,10);
metrics = [0.068, 0, 0];
usemetrics = [true false false]; % only fit the distal density

fitkye = zeros(size(kpylist2));
YPrat2 = fitkye;
fitkyeparams = cell(size(kpylist2));
fitdens = fitkyeparams;

for kc = 1:length(kpylist2)
  params.kpy = kpylist2(kc); % lysosome production rate
  
  fitguess = 1;
  minopt = optimset('Display','final','TolFun',1e-4,'TolX',1e-4,'UseParallel',true);
  lb = 0.01; % lower bounds
  ub = 100; % upper bounds
  fitkye(kc) = fmincon(@(fitp) fitfunc_tipfusion([fitparams(1) params.kpy fitp],metrics,usemetrics,params,options), fitguess,[],[],[],[],lb,ub,[],minopt);
  
  [res, gotmetrics,fitkyeparams{kc},fitdens{kc}] = fitfunc_tipfusion([params.pfra kpylist2(kc) fitkye(kc)],metrics,[true true true],params,options);
  YPrat2(kc) = gotmetrics(3);
  disp([kc,YPrat2(kc)])
end
save('./data/onedone_nobranch_fitkpykye.mat')
%% plot proximal density ratio vs kpy
scl = params.vpreal/params.Lreal*60;
YPratexp = 1.88;
% figure 6 panel b
figure('Name','Figure 6(b)','NumberTitle','off')
plot(kpylist2*scl,YPrat2,'b','LineWidth',2)
hold all
tmp = [min(kpylist2),max(kpylist2)]*scl;
xlim(tmp)
plot(tmp,YPratexp*[1,1],'k--')
hold off
plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28)
title(sprintf('fusion probability p_f = %0.03f',fitparams(1)))
xlabel('lysosome production rate k_p^y (min^{-1})')
ylabel('ratio of proximal LAMP1/LC3 density')
%% compute densities for fitted parameters
params.pf = fitparams(1);
% first letter denotes AV direction, second letter denotes lysosome direction
% a: anterograde, r: retrograde, s: stationary
% fusion probability is currently assumed to be independent of the
% direction of relative motion
params.pfra = params.pf; params.pfrr = params.pf; params.pfrs = params.pf;
params.pfaa = params.pf; params.pfar = params.pf; params.pfas = params.pf;
params.pfsa = params.pf; params.pfsr = params.pf;
params.pft = params.pf; % fusion probability within tip

params.kpy = fitparams(2);
params.kye = fitparams(3);

[dens,success] = onedone_branches_solver(params,options);

% plots for densities with fully fitted values
xvals = dens.xvals;
dx = xvals(2)-xvals(1);

% plot distributions
Btot = dens.Ba+dens.Bs+dens.Br;
Ptot = dens.S + dens.R +Btot;
Stot = Btot;
Ptf0 = dens.Pt0 - dens.Ptu0;

Ytot = dens.Ya+dens.Yr+dens.Ysa+dens.Ysr+dens.Yuw+dens.Yus;
Butot = dens.Bau+dens.Bsu+dens.Bru;
Putot = dens.Su+dens.Ru+Butot;
Ftot = Ptot-Putot;
YFtot = Ytot+Ftot; % LAMP1 puncta;
%% plot distribution of fused AVs
% figure 6 panel c
figure('Name','Figure 6(c)','NumberTitle','off')
plot(xvals*params.Lreal,Ftot./Ptot,'b','LineWidth',2)
plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28);
xlabel('position (\mum)')
ylabel('fraction of AVs fused')
ylim([0,1])
xlim([0,1060])
%% panels d,e, f require running simulations using code provided at 
%% https://github.com/lenafabr/particleDynamics1D
%% read data for panel d
fnamestr = '../results/snapfiles/getmoviesnaps_nobranch_onedone.snap.out';
options = struct();
options.getmoviesnaps = 1;
[grouplist,tvals] = readsnapshot(fnamestr,options);
save('./data/getSimSnapshot_nobranch_onedone.mat')
%% alternatively, load the mat file directly
load('./data/getSimSnapshot_nobranch_onedone.mat')
%% plot a snapshot from simulations
snaptime = 5; % time at which snapshot is taken (max 10) in non-dimensionalized units (Lreal/vpreal)
% snaptime should be greater than ~1 for steady state
figure('Name','Figure 6(d)','NumberTitle','off')
plotSimSnapShot(snaptime,grouplist,tvals)
%% read data from simulation snap file
fnamestr = './data/matchexptparams_nobranch_onedone.snap.out';
options = struct('readlastsnap',1);

[grouplist,tvals,domlen,ntrials] = readsnapshot(fnamestr,options);

% read phagosome positions
phagopos = [];
lysopos = [];
phxfuse = [];

for tc = 1:length(grouplist)
	pos = vertcat(grouplist(tc).particles(:).pos);
	ptypes = vertcat(grouplist(tc).particles(:).type);
	xfuse = vertcat(grouplist(tc).particles(:).xfuse);
  intip = vertcat(grouplist(tc).particles(:).intip);
	
	phagopos = cat(1,phagopos,pos(ptypes==2));
	lysopos = cat(1,lysopos,pos(ptypes==1 & intip==0));
	phxfuse = cat(1,phxfuse,xfuse(ptypes==2));
end
save('./data/matchexptparams_nobranch_onedone.mat');
%% alternatively, load the mat file directly
load('./data/matchexptparams_nobranch_onedone.mat');
%% plot organelle densities
vp = 1;
vy = 2.67;
Lreal = params.Lreal;
dt = 1d-4; % simulation timestep
phfrate = 350.001; % framerate for phagosomes (binwidth)
lyfrate = 100.001; % framerate for lysosomes (binwidth)
xp = 0:(phfrate*vp*dt):1; % possible positions for phagosomes
xp(end) = 1;
xl = 0:(lyfrate*vy*dt):1; % possible positions for lysosomes
xl(end) = 1;
% organelle densities from sims
[nphago,~,bininds] = histcounts(phagopos,xp);
nfusehist = accumarray(bininds,phnfuse)'./nphago;
xpvals = 0.5*(xp(1:end-1)+xp(2:end));
dxp = (xpvals(2)-xpvals(1));

cphago = [215,38,156]/255;
clyso = [126,212,238]/255;
% figure 6 panel e
figure('Name','Figure 6(e)','NumberTitle','off')
plot(xpvals*Lreal,fliplr(nphago)./ntrials./dxp./Lreal,'o','Color',cphago,'MarkerFaceColor',cphago)
nfusedphago = histcounts(phagopos(phnfuse>0),xl);
nlyso = histcounts(lysopos,xl);
xlvals = 0.5*(xl(1:end-1)+xl(2:end));
dxl = (xlvals(2)-xlvals(1));

hold on
plot(xlvals*Lreal,fliplr(nlyso+nfusedphago)./ntrials./dxl./Lreal,'o','Color',clyso,'MarkerFaceColor',clyso)
plot(xvals*Lreal,Ptot/Lreal,'Color',cphago)
plot(xvals*Lreal,YFtot/Lreal,'Color',clyso)
hold off

legend('AV density, simulations','Lysosome density, simulations',...
			'AV density, analytical model','Lysosome density, analytical model')

plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28)
xlim([0,1050]);
ylabel('organelle density (\mum^{-1})')
xlabel('distance from distal end (\mum)')
%% plot histogram of position of first fusion
% figure 6 panel f
figure('Name','Figure 6(f)','NumberTitle','off')
binwidth = 10;
% xfusevals = (1-phxfuse(phxfuse>0 & phxfuse<1))*params.Lreal; % fused in the domain
xfusevals = (1-phxfuse(phxfuse>0))*params.Lreal; % fused anywhere
h = histogram(xfusevals,'BinWidth',binwidth,'Normalization','probability');
xlim([0,100])
xlabel('distance from distal tip at the time of fusion (\mum)')
ylabel('fraction of fused AVs')
plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28)