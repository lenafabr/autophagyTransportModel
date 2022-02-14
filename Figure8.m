%% script to reproduce plots in figure 6
% requires running Figure7.m to generate a savefile with fitted parameters
% required files: fitparams_onedone_branches.mat
clear;clc;
%% add scripts to path
addpath('./tools/')
%% set up options for calculation
options = struct();

options.nx = 1000; % number of points on the calculation grid

options.solvePfuse = true;
options.dotipfusion = true; % fuse particles in the distal tip
options.onedone = true; % allow AVs to only fuse once

numBranch = 5; % number of branches

xvals = linspace(0,1,options.nx);
dx = xvals(2)-xvals(1);
%% load fitted parameters
load('./data/fitparams_onedone_branches.mat','fitparams','params','dens')
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

xvals = dens.xvals;
dx = xvals(2)-xvals(1);

Btot = dens.Ba+dens.Bs+dens.Br;
Ptot = dens.S + dens.R +Btot;
Stot = Btot;
Ptf0 = dens.Pt0 - dens.Ptu0;

Ytot = dens.Ya+dens.Yr+dens.Ysa+dens.Ysr+dens.Yuw+dens.Yus;
Butot = dens.Bau+dens.Bsu+dens.Bru;
Putot = dens.Su+dens.Ru+Butot;
Ftot = Ptot-Putot;
YFtot = Ytot+Ftot; % LAMP1 puncta;
%% Scan over different kd
kdlistreal = logspace(-5,-1.2,30); % per min
kdlist = kdlistreal/params.vpreal*params.Lreal;

xdist = 249/params.Lreal;
xmid = 251/params.Lreal;

distind = find(xvals<xmid);
midind = find(xvals>xmid & xvals<1-xmid);
proxind = find(xvals>(1-xdist));

Ptotdist = (sum(Ptot(distind))-0.5*(Ptot(distind(1))+Ptot(distind(end))))*dx;
Pprox = (sum(Ptot(proxind))-0.5*(Ptot(proxind(1))+Ptot(proxind(end))))*dx;
Pmid = (sum(Ptot(midind))-0.5*(Ptot(midind(1))+Ptot(midind(end))))*dx;
clear fracdegdist fracdegmid fracdegprox
for kc = 1:length(kdlist)
  params.kd = kdlist(kc);
  densdeg = calcDecayedAVdistribution(params,dens);
  alldensdeg(kc) = densdeg;
  
  totundeg = densdeg.Rf+densdeg.Brf+densdeg.Baf+densdeg.Bsf+densdeg.Sf;
  totdeg = Ftot - totundeg;
  
  degdist = (sum(totdeg(distind))-0.5*(totdeg(distind(1))+totdeg(distind(end))))*dx;
  degmid = (sum(totdeg(midind))-0.5*(totdeg(midind(1))+totdeg(midind(end))))*dx;
  degprox = (sum(totdeg(proxind))-0.5*(totdeg(proxind(1))+totdeg(proxind(end))))*dx;
  
  fracdegdist(kc) = degdist./Ptotdist;
  fracdegprox(kc) = degprox./Pprox;
  fracdegmid(kc) = degmid/Pmid;
end

%% make plots
cmap = BBVYWcolormap(5);
cmap = cmap(2:end-1,:);
taulist = 1./kdlist/params.vpreal*params.Lreal/60;
distind = find(fracdegdist<1 & fracdegdist>0);
proxind = find(fracdegprox<1 & fracdegprox > 0);
midind = find(fracdegmid<1 & fracdegmid>0);
% figure 8 panel g
figure('Name','Figure 8(g)','NumberTitle','off')
semilogx(taulist(distind),fracdegdist(distind),'Color',cmap(1,:))
hold on
semilogx(taulist(midind),fracdegmid(midind),'Color',cmap(2,:))
semilogx(taulist(proxind),fracdegprox(proxind),'Color',cmap(3,:))
ylim([0,1])
xlim([0.4,max(taulist)])

plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28)
xlabel('degradation timescale 1/k_d (min)')
ylabel('fraction of AVs with fully degraded IAM')

% experimental measurements
hold all
fdhippo = [0.084 0.274 0.19];
stdhippo = [0.033 0.065 0.06];
thippo(1) = interp1(fracdegdist,taulist,fdhippo(1));
xrangehippo(1,:) = interp1(fracdegdist,taulist,fdhippo(1)+stdhippo(1)*[1,-1]);
thippo(2) = interp1(fracdegprox,taulist,fdhippo(2));
xrangehippo(2,:) = interp1(fracdegprox,taulist,fdhippo(2)+stdhippo(2)*[1,-1]);
thippo(3) = interp1(fracdegmid,taulist,fdhippo(3));
xrangehippo(3,:) = interp1(fracdegmid,taulist,fdhippo(3)+stdhippo(3)*[1,-1]);
hippocol = [0.8,0.7,0];
h1 = errorbar(thippo,fdhippo,stdhippo,stdhippo,thippo-xrangehippo(:,1)', xrangehippo(:,2)'-thippo,...
  '.','MarkerSize',27,'Color',hippocol,'LineWidth',2);

hold all
stddrg = [0.0375 0.055];
fddrg = [0.4825 0.69];
tdrg(1) = interp1(fracdegdist,taulist,fddrg(1));
xrangedrg(1,:) = interp1(fracdegdist,taulist,fddrg(1)+stddrg(1)*[1,-1]);
tdrg(2) = interp1(fracdegprox,taulist,fddrg(2));
xrangedrg(2,:) = interp1(fracdegprox,taulist,fddrg(2)+stddrg(2)*[1,-1]);
drgcol = [0,0.7,0];
h2 = errorbar(tdrg,fddrg,stddrg,stddrg,tdrg-xrangedrg(:,1)', xrangedrg(:,2)'-tdrg,...
  '.','MarkerSize',27,'Color',drgcol,'LineWidth',2);
hold off
legend('distal','mid','proximal','measured: hippocampal','measured: DRG','FontSize',27)
ax = gca;
ax.XTick = [1,10,100,1000];
%% plot distribution of degraded AVs
tmean = mean(thippo);
kdmean = 1/(tmean*60)*params.Lreal/params.vpreal;
params.kd = kdmean;
densdeg = calcDecayedAVdistribution(params,dens);

% total fused organelles
Ftot = dens.R+dens.Br+dens.Ba+dens.Bs+dens.S - dens.Bau-dens.Bru-dens.Ru-dens.Bsu-dens.Su;

% plot fraction degraded
totundeg = densdeg.Rf+densdeg.Brf+densdeg.Baf+densdeg.Bsf+densdeg.Sf;
totdeg = Ftot - totundeg;
fracdeg = totdeg./Ptot;
% figure 8 panel (h)
figure('Name','Figure 8(g)','NumberTitle','off')
plot(xvals*params.Lreal,fracdeg,'.-')
xlim([0,1055])
ax = gca;
ax.YTick = [0:0.05:0.3];
plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28);
xlabel('position (\mum)')
ylabel('fraction of AVs with degraded IAM')