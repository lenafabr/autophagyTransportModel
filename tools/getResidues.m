function [res, gotmetrics,params,dens] = getResidues(fitparams,metrics,usemetrics,params,options)
% fitting function gives residuals to desired metrics
% fit params: [pf, kpy, kye]
% metrics: [rhoy, distal frac fused, proximal Y:P ratio]
% usemetrics = logical for which metrics are used for final residual

opt.onedone = false;
opt.scaleres = false;

if (exist('options','var'))
    opt = copyStruct(options,opt,'addnew',true);
end

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

if(opt.onedone)
  [dens,success] = onedone_branches_solver(params,opt);
else
  [dens,success] = unlim_branches_solver(params,opt);
end

if ~success
    res = NaN;
    return
end

xvals = dens.xvals;
dx = xvals(2)-xvals(1);
xd = opt.xd;

% plot distributions
Btot = dens.Ba+dens.Bs+dens.Br;
Ptot = dens.R+Btot+dens.S;

Ytot = dens.Ya+dens.Yr+dens.Ysa+dens.Ysr+dens.Yuw+dens.Yus;

Butot = dens.Bau+dens.Bsu+dens.Bru;
Putot = dens.Ru+Butot+dens.Su;
Ftot = Ptot-Putot;
YFtot = Ytot+Ftot; % LAMP1 puncta;

% distal fraction fused
distind = find(xvals<xd);

Ptotdist = (sum(Ptot(distind))-0.5*(Ptot(distind(1))+Ptot(distind(end))))*dx;
Fdist =  (sum(Ftot(distind))-0.5*(Ftot(distind(1))+Ftot(distind(end))))*dx;
Ffdist = (Fdist)./(Ptotdist);

% distal LC3 density
Pdist = (Ptotdist)/xd/params.Lreal;
% distal LAMP1 density
Ytotdist = (sum(Ytot(distind))-0.5*(Ytot(distind(1))+Ytot(distind(end))))*dx;
YFtotdist = (sum(YFtot(distind))-0.5*(YFtot(distind(1))+YFtot(distind(end))))*dx;
Ydist = (YFtotdist+dens.Yt0)/xd/params.Lreal;

proxind = find(xvals>(1-xd));
YPrat = sum(YFtot(proxind))/sum(Ptot(proxind));

gotmetrics = [Ydist,Ffdist,YPrat];

if(opt.scaleres)
  res = norm((metrics(usemetrics) - gotmetrics(usemetrics))./metrics(usemetrics));
else
  res = norm(metrics(usemetrics) - gotmetrics(usemetrics));
end

end