function [dens,success,params,sol] = onedone_branches_solver(params,options)
% obtains AV and lysosome densities for the one and done fusion model
% all lengths are scaled by the length of the domain: params.Lreal
% all times are scaled by time required by an AV to cross the domain: params.Lreal/params.vpreal

opt = struct();
opt.nx = 1000;
opt.dotipfusion = false;

if (exist('options','var'))
  opt = copyStruct(options,opt,'addnew',1);
end

nb = length(params.branchpoints);

%% set up fusion rates in params structure
phit = ones(1,nb+1); phitp = ones(1,nb+1);
for c = nb+1: -1 : 1
  
  % AV retrograde, lyso anterograde
  tpass = params.ell/(params.vp+params.vya);
  params.kfra(c) = params.pfra*phit(c)/tpass;    params.kfrap(c) = params.pfra*phitp(c)/tpass;
  tpass = params.ell/(params.vpb+params.vya);
  params.kfrab(c) = params.pfra*phit(c)/tpass;   params.kfrapb(c) = params.pfra*phitp(c)/tpass;
  
  % AV anterograde, lyso retrograde
  tpass = params.ell/(params.vp+params.vyr);
  params.kfar(c) = params.pfar*phit(c)/tpass;    params.kfarp(c) = params.pfar*phitp(c)/tpass;
  tpass = params.ell/(params.vpb+params.vyr);
  params.kfarb(c) = params.pfar*phit(c)/tpass;   params.kfarpb(c) = params.pfar*phitp(c)/tpass;
  
  % AV retrograde, lyso halted
  tpass = params.ell/(params.vp);
  params.kfrs(c) = params.pfrs*phit(c)/tpass;    params.kfrsp(c) = params.pfrs*phitp(c)/tpass;
  % AV anterograde, lyso halted
  params.kfas(c) = params.pfas*phit(c)/tpass;    params.kfasp(c) = params.pfas*phitp(c)/tpass;
  
  tpass = params.ell/(params.vpb);
  params.kfrsb(c) = params.pfrs*phit(c)/tpass;   params.kfrspb(c) = params.pfrs*phitp(c)/tpass;
  params.kfasb(c) = params.pfas*phit(c)/tpass;   params.kfaspb(c) = params.pfas*phitp(c)/tpass;
  
  
  % moving in same direction retrograde
  tpass = params.ell/abs(params.vyr - params.vp);
  params.kfrr(c) = params.pfrr*phit(c)/tpass;    params.kfrrp(c) = params.pfrr*phitp(c)/tpass;
  tpass = params.ell/abs(params.vyr - params.vpb);
  params.kfrrb(c) = params.pfrr*phit(c)/tpass;   params.kfrrpb(c) = params.pfrr*phitp(c)/tpass;
  
  % moving in same direction anterograde
  tpass = params.ell/abs(params.vya - params.vp);
  params.kfaa(c) = params.pfaa*phit(c)/tpass;    params.kfaap(c) = params.pfaa*phitp(c)/tpass;
  tpass = params.ell/abs(params.vya - params.vpb);
  params.kfaab(c) = params.pfaa*phit(c)/tpass;   params.kfaapb(c) = params.pfaa*phitp(c)/tpass;
  
  % AV stopped, lyso anterograde
  tpass = params.ell/(params.vya);
  params.kfsa(c) = params.pfsa*phit(c)/tpass;    params.kfsap(c) = params.pfsa*phitp(c)/tpass;
  
  % AV stopped, lyso retrograde
  tpass = params.ell/(params.vyr);
  params.kfsr(c) = params.pfsr*phit(c)/tpass;    params.kfsrp(c) = params.pfsr*phitp(c)/tpass;
end

%% solve for total phagosome densities (no fusion involved, can handle branches)
dens = calcTotalAVdistribution(params,opt);
%% solve for unfused phagosomes and lysosomes using matlab's built-in PDE solver
xmesh = linspace(0,1,round(length(dens.xvals)/(nb+1)));
solinit = bvpinit(xmesh, @(x) init_guess_onedone(params));
%% run solver
soloptions = bvpset('NMax',1000);
sol = bvp4c(@(x,y) diffeqs_onedone(y,params), @(ya,yb) boundcond_onedone(ya,yb,params), solinit,soloptions);

%% get tip lysosomes
Yt0 = params.vya*sol.y(4,1)/(params.kye + params.pft*(params.vpb*sol.y(2,1) + params.kpp0));

Yt = zeros(1,nb);
for c = 1:nb
  sc = 5*(nb+c);
  Yt(c) = params.vya*sol.y(sc+4,1)/(params.kye + params.pft*(params.vpb*sol.y(sc+2,1) + params.kpp(c)));
end
%% rescale solution
branchlen = diff([0 params.branchpoints params.L]);
branchlenp = params.branchlen;
branchpoints = [0 params.branchpoints params.L];

for c = 1:nb+1
  ind = dens.mainind{c};
  xshift = dens.xvals(ind) - branchpoints(c);
  
  sc = 5*(c-1);
  dens.Bru(ind) = interp1(sol.x*branchlen(c),sol.y(sc+1,:),xshift);
  dens.Bau(ind) = interp1(sol.x*branchlen(c),sol.y(sc+2,:),xshift);
  
  dens.Ru(ind) = interp1(sol.x*branchlen(c),sol.y(sc+3,:),xshift);
  dens.Ya(ind) = interp1(sol.x*branchlen(c),sol.y(sc+4,:),xshift);
  dens.Yr(ind) = interp1(sol.x*branchlen(c),sol.y(sc+5,:),xshift);
  
  
  fs = params.ell*(params.kfsa(c)*dens.Ya(ind) + params.kfsr(c)*dens.Yr(ind));
  dens.Bsu(ind) = params.khb*(dens.Bau(ind)+dens.Bru(ind))./(params.kwb + params.ks + fs);
  dens.Su(ind) = params.khr*(dens.Ru(ind))./(params.kwr + fs);
  
  fys = params.ell*(params.kfasb(c)*dens.Bau(ind) + params.kfrsb(c)*dens.Bru(ind) + params.kfrs(c)*dens.Ru(ind));
  dens.Ysa(ind) = params.khy(1)*dens.Ya(ind)/(params.kwy + fys);
  dens.Ysr(ind) = params.khy(1)*dens.Yr(ind)/(params.kwy + fys);
end

for c = 1:nb
  sc = 5*(nb+c);
  ind = dens.collind{c};
  xcoll = dens.xvals(ind);
  dens.Brup{c} = interp1(sol.x*branchlenp(c),sol.y(sc+1,:),xcoll);
  dens.Baup{c} = interp1(sol.x*branchlenp(c),sol.y(sc+2,:),xcoll);
  dens.Rup{c} = interp1(sol.x*branchlenp(c),sol.y(sc+3,:),xcoll);
  dens.Yap{c} = interp1(sol.x*branchlenp(c),sol.y(sc+4,:),xcoll);
  dens.Yrp{c} = interp1(sol.x*branchlenp(c),sol.y(sc+5,:),xcoll);
  
  fs = params.ell*(params.kfsap(c)*dens.Yap{c} + params.kfsrp(c)*dens.Yrp{c});
  dens.Bsup{c} = params.khb*(dens.Baup{c}+dens.Brup{c})./(params.kwb + params.ks + fs);
  dens.Sup{c} = params.khr*(dens.Rup{c})./(params.kwr + fs);
  
  fys = params.ell*(params.kfaspb(c)*dens.Baup{c} + params.kfrspb(c)*dens.Brup{c} + params.kfrsp(c)*dens.Rup{c});
  dens.Ysap{c} = params.khy(1)*dens.Yap{c}/(params.kwy + fys);
  dens.Ysrp{c} = params.khy(1)*dens.Yrp{c}/(params.kwy + fys);
end

dens.kfra = params.kfra;
dens.kfar = params.kfar;
dens.kfrs = params.kfrs;
dens.kfas = params.kfas;
dens.kfrr = params.kfrr;
dens.kfaa = params.kfaa;
dens.kfsa = params.kfsa;
dens.kfsr = params.kfsr;

dens.kfrap = params.kfrap;
dens.kfarp = params.kfarp;
dens.kfrsp = params.kfrsp;
dens.kfasp = params.kfasp;
dens.kfrrp = params.kfrrp;
dens.kfaap = params.kfaap;
dens.kfsap = params.kfsap;
dens.kfsrp = params.kfsrp;

dens.phit = phit;

dens.Yt0 = Yt0;
dens.Yt = Yt;
dens.Pt0 = 0;
dens.Ptu0 = 0;

% unused fields
dens.Yuw = 0*dens.Bru;
dens.Yus = dens.Yuw;

success = true;

end