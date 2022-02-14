%% script to reproduce plots in figure 7 Supplement 2
clear;clc;
%% add scripts to path
addpath('./tools/')
%% set up options for calculation
options = struct();

options.nx = 1000; % number of points on the calculation grid

options.solvePfuse = true;
options.dotipfusion = true; % fuse particles in the distal tip
options.onedone = false; % allow AVs to only fuse once

numBranch = 5; % number of branches

xvals = linspace(0,1,options.nx);
dx = xvals(2)-xvals(1);
%% set up parameter structure
params = struct();

% domain length and time scales
params.Lreal = 1055; % length of domain in real units (um). Sets length scale
params.vpreal = 0.75; % AV velocity in real units (um/s). Sets time scale along with Lreal
params.L = 1; % scaled domain length
params.vp = 1; % scaled AV velocity

% organelle sizes
params.rp = 0.4/params.Lreal; % AV radius
params.ry = 0.1/params.Lreal; % lysosome radius
params.ell = 2*(params.rp+params.ry); % minimum contact distance

% AV production params
params.kpp0 = 6.97; % AV production rate in main axon tip
params.kpp = params.kpp0*ones(1,numBranch); % AV production rate in branch tips

% params for bidirectional AV motion
params.vpb = 1; % AV velocity in bidirectional state
lam = 2; % AV run length in bidirectional state (um)
run_length = lam/params.Lreal; % scaled run length
fstat = 0.4603; % fraction of stationary AVs (obtained from exp data)
params.khb = params.vpb/run_length; % rate of halting in the bidirectional state
params.kwb = params.khb*(1-fstat)/fstat; % rate of walking in the bidirectional state
params.pretro = 0.5; % probability of walking retrograde upon restarting

% params for persistently retrograde AV motion
params.ks = 1; % rate of AVs switching to persistent retrograde motion
params.kwr = 1; % rate of walking in the retrograde state
haltratio = 0.5967; % ratio of halting rate to walking rate for retrograde AVs
params.khr = haltratio*params.kwr; % rate of halting in the retrograde state

% parameters for lysosome production and motion
params.kpy = 78.92; % lysosome production rate
% 40% retro @ 1.44um/s, 30% stationary, 30% antero @ 2um/s
% source: Boecker et al. 2020
params.vya = 0.3*2/(0.3+0.15)/params.vpreal; % anterograde lysosome velocity
params.vyr = 0.4*1.5/(0.4+0.15)/params.vpreal; % retrograde lysosome velocity

params.khy = 0; %lysosome halting rate
params.kwy = 0; %lysosome walking rate

% parameters for distal tip
params.kye = 0.97; % rate of lysosome return from distal tip
params.kpe = inf; % rate of AV return from distal tip

% parameters for lysosome decay
params.kda = 0; params.kdr = 0; params.kds = 0;

% parameters for fusion probability
params.pf = 0.0074;
% first letter denotes AV direction, second letter denotes lysosome direction
% a: anterograde, r: retrograde, s: stationary
% fusion probability is currently assumed to be independent of the
% direction of relative motion
params.pfra = params.pf; params.pfrr = params.pf; params.pfrs = params.pf;
params.pfaa = params.pf; params.pfar = params.pf; params.pfas = params.pf;
params.pfsa = params.pf; params.pfsr = params.pf;
params.pft = params.pf; % fusion probability within tip

% options and parameters for axonal branching
options.pfsplit = false; %same fusion probability in all branches
params.mainrad = 1*ones(1,numBranch+1); % radii of main axon segments
params.collrad = 1*ones(1,numBranch+1); % radii of axon collaterals/branches
params.branchpoints = linspace(250,params.Lreal-250,numBranch)/params.Lreal; % location of branch points along the main axon
params.branchlen = 164/params.Lreal*ones(size(params.branchpoints)); % branch lengths

% split organelles according to number of branch tips served
params.psplit = [];
for c = 1:numBranch
  params.psplit(c) = c/(c+1);
end

%% fraction of retrograde AVs vs switching timescale
kpp0save = params.kpp0;
tswitch = linspace(0.01,60,100); % switching timescale in minutes
kslist = 1./tswitch/60/params.vpreal*params.Lreal; % switching rate

Frdist = zeros(size(kslist)); % distal fraction retrograde
Frmid = Frdist; % mid axon fraction retrograde
Pdist = Frdist; % distal AV density

xdist = 249/params.Lreal;
xmid = 251/params.Lreal;

for kc = 1:length(kslist)
  params.ks = kslist(kc);
  
  dens = calcTotalAVdistribution(params,options);
  Btot = dens.Ba+dens.Bs+dens.Br;
  Ptot = dens.R+Btot+dens.S;
  
  % distal fraction retrograde
  distind = find(xvals<xdist);
  Ptotdist = (sum(Ptot(distind))-0.5*(Ptot(distind(1))+Ptot(distind(end))))*dx;
  Rdist = (sum(dens.R(distind))-0.5*(dens.R(distind(1))+dens.R(distind(end))))*dx;
  Frdist(kc) = Rdist/Ptotdist;
  
  % mid fraction retro
  midind = find(xvals>xmid & xvals < 1-xmid);
  Rmid = (sum(dens.R(midind))-0.5*(dens.R(midind(1))+dens.R(midind(end))))*dx;
  Pmid = (sum(Ptot(midind))-0.5*(Ptot(midind(1))+Ptot(midind(end))))*dx;
  Frmid(kc) = Rmid/Pmid;
  
  % distal density
  Pdist(kc) = Ptotdist/xdist/params.Lreal;
end

Frexpdist = 0.124;
Frexpmid = 0.602;

ind = find(~isnan(Frdist));
ksfit = interp1(Frdist(ind),kslist(ind),Frexpdist);
ksmid = interp1(Frmid(ind),kslist(ind),Frexpmid);
err = abs(ksfit-ksmid)/ksfit;
params.ks = ksfit;
fprintf('Fitted switching rate k_s = %0.2g\n',ksfit)
%% What production rate is necessary to give correct distal density?
Pdistexp = 0.045; % per um
Pdistfit = interp1(kslist,Pdist,ksfit);
kppfit = Pdistexp/Pdistfit*kpp0save;
params.kpp = kppfit*ones(1,numBranch);
params.kpp0 = kppfit;
fprintf('Fitted AV production rate k_p^p = %0.2g\n',kppfit)
%% directly fit 3 parameters : pf, kpy, kye
% vs 3 metrics: distal lysosome density, distal fraction fused, proximal Y:P ratio
metrics = [0.068 0.49 1.88];
usemetrics = true(1,3);  % which metrics to use in calculations

% options for fitting
options.xd = xdist;
options.nx = 1000;
options.scaleres = true; 

% for model with no branches
fitguess = [[0.008,104,0.72]];
[res,gotmetrics] = getResidues(fitguess,metrics,usemetrics,params,options);
fprintf('3 parameter fitting\n')
fprintf('residue: %0.3g\n',res)
fprintf('metrics: %0.2g, %0.2g, %0.3g\n', gotmetrics)
%% run matlab minimizer
% change 'PlotFcn' to 'optimplotfval' to visualize the fitting process
minopt = optimset('Display','none','TolFun',1e-6,'TolX',1e-6,...
  'UseParallel',true,'PlotFcn','optimplotfval');
lb = [0,0,0]; % lower bounds
ub = [1,Inf,Inf]; % upper bounds
fprintf('Running 3 parameter fit...\n')
fitparams = fmincon(@(fitp) getResidues(fitp,metrics,usemetrics,params,options), fitguess,[],[],[],[],lb,ub,[],minopt);
%% check fitted values
[res, gotmetrics,params,dens] = getResidues(fitparams,metrics,usemetrics,params,options);
fprintf('Obtained metrics:\n')
fprintf('distal lyso density: %0.3f\n',gotmetrics(1))
fprintf('distal AV fraction fused: %0.2f\n',gotmetrics(2))
fprintf('proximal Y:P ratio: %0.2f\n',gotmetrics(3))
%% save fitted parameters
save('./data/fitparams_unlim_branches.mat','fitparams','params','dens')
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
    
    [dens,success] = onedeon_branches_solver(params,options);
    
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
  save('./data/unlim_branches_scanpfkpy_kyefit.mat')
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
% figure 7S2 panel a
figure('Name','Figure 7 Supplement 2(a)','NumberTitle','off')
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
save('./data/unlim_branches_fitkpykye.mat')
%% plot proximal density ratio vs kpy
scl = params.vpreal/params.Lreal*60;
YPratexp = 1.88;
% figure 7S2 panel b
figure('Name','Figure 7 Supplement 2(b)','NumberTitle','off')
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

[dens,success] = unlim_branches_solver(params,options);

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
%% average analytical distribution over bins
bw = 10;
nint = round(bw/params.Lreal/dx); % integrate over this many points
clear xavg Ptotavg Btotavg LAMPavg
LAMP = Ytot+Ftot;
for bc = 1:floor(length(xvals)/nint)-1
  ind = nint*(bc-1)+1;
  xavg(bc) = mean(xvals(ind:ind+nint));
  Ptotavg(bc) = mean(Ptot(ind:ind+nint));
  LAMPavg(bc) = mean(LAMP(ind:ind+nint));
end
% figure 7 supplement 2 panel c
figure('Name','Figure 7 Supplement 2(c)','NumberTitle','off')
plot(xavg*params.Lreal,Ptotavg/params.Lreal,'m-.','LineWidth',2)
hold all
plot(xavg*params.Lreal,LAMPavg/params.Lreal,'-.','LineWidth',2,'Color','c')
hold off
legend('LC3 density','LAMP1 density','FontSize',28,'Location','northeast')
xlim([0,100])
xlabel('position (\mum)')
ylabel('density (\mum^{-1})')
plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28);
%% plot number of fusions and distribution of fused AVs
fd = numfusion_branches(params,dens);
avgNf = (fd.Baf + fd.Bsf+fd.Brf + fd.Rf+fd.Sf)./(dens.Ba + dens.Bs + dens.Br + dens.R + dens.S);

% figure 7 supplement 2 panel d
figure('Name','Figure 7 Supplement 2(d)','NumberTitle','off')

xlabel('position (\mum)')
xlim([0,1050])

yyaxis left
plot(xvals*params.Lreal,Ftot./Ptot,'b','LineWidth',2)
plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28);
ylabel('fraction of AVs fused')
ylim([0,1])

yyaxis right
plot(xvals*params.Lreal,avgNf,'r','LineWidth',2)
plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28);
ylabel('average number of fusions per AV')