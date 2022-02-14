function fd = calcDecayedAVdistribution(params,d)
% given output densities from solver
% calculate the density of fusions (ie: # fusion * AV density) in main axon
% and branches
% d = intput densities structure
% returns structure containing all densities

xvals = d.xvals;
dx = xvals(2)-xvals(1);
numBranch = length(params.branchpoints);
branchpoints = [0 params.branchpoints params.L];
%% on main axon, bidirectional

Mf(1,1) = -(params.ks + params.kd + params.khb)+ params.khb*params.kwb*params.pretro/(params.ks + params.kd + params.kwb);
Mf(1,2) = params.khb*params.kwb*params.pretro/(params.ks + params.kd+params.kwb);
Mf(2,1) = - params.khb*params.kwb*(1-params.pretro)/(params.ks + params.kd + params.kwb);
Mf(2,2) = (params.ks + params.kd+params.khb)- params.khb*params.kwb*(1-params.pretro)/(params.ks + params.kd + params.kwb);
Mf = Mf/params.vp;

[evecf,evalf] = eig(Mf);
evalf = diag(evalf);

%%
mainind = {};
for c = 1:numBranch+1
  ind = find(xvals>=branchpoints(c) & xvals<=branchpoints(c+1));
  mainind{c} = ind;
  xshift = xvals(ind) - branchpoints(c);
  
  fbr = params.ell*d.Bru(ind).*(d.kfra(c).*d.Ya(ind) + d.kfrs(c).*(d.Ysa(ind)+d.Ysr(ind)));
  fbs = params.ell*d.Bsu(ind).*(d.kfsa(c).*d.Ya(ind) + d.kfsr(c).*d.Yr(ind));
  fba = params.ell*d.Bau(ind).*(d.kfar(c).*d.Yr(ind) + d.kfas(c).*(d.Ysa(ind)+d.Ysr(ind)));
  
  gvals = [fbr + params.kwb*params.pretro/(params.ks + params.kd + params.kwb)*fbs; ...
    -fba - params.kwb*(1-params.pretro)/(params.ks+ params.kd + params.kwb)*fbs];
  gvals = gvals/params.vp;
  
  
  Vinvg = zeros(2,length(ind));
  for xc = 1:length(xshift)
    Vmat = [evecf(:,1)*exp(evalf(1)*xshift(xc)) , evecf(:,2)*exp(evalf(2)*xshift(xc))];
    Vinvg(:,xc) = Vmat\gvals(:,xc);
  end
  
  for vc = 1:2
    zvals(vc,ind) = (cumsum(Vinvg(vc,:)) - 0.5*(Vinvg(vc,1)+Vinvg(vc,:)))*dx;
  end
  
  dzvals(:,c) = zvals(:,ind(end))-zvals(:,ind(1));
end


%% along collateral branches, solve for bidirectional densities
% v1*exp(l1*x)*z1 + v2*exp(l2*x)*z2 evaluated at end of collateral

zvalsp = {}; vzvalsp={};
for c = 1:numBranch
  %%
  collind{c} = find(xvals<=params.branchlen(c));
  ind = collind{c};
  xcoll = xvals(ind);
  
  fbr = params.ell*d.Brup{c}.*(d.kfrap(c).*d.Yap{c} + d.kfrsp(c).*(d.Ysap{c}+d.Ysrp{c}));
  fbs = params.ell*d.Bsup{c}.*(d.kfsap(c).*d.Yap{c} + d.kfsrp(c).*d.Yrp{c});
  fba = params.ell*d.Baup{c}.*(d.kfarp(c).*d.Yrp{c} + d.kfasp(c).*(d.Ysap{c}+d.Ysrp{c}));
  
  %%
  gvalsp{c} = [fbr + params.kwb*params.pretro/(params.ks + params.kd+params.kwb)*fbs; ...
    -fba - params.kwb*(1-params.pretro)/(params.ks+params.kd+params.kwb)*fbs];
  gvalsp{c} = gvalsp{c}/params.vp;
  
  Vinvgp{c} = zeros(size(gvalsp{c}));
  for xc = 1:length(xcoll)
    xcoll = xvals(ind);
    Vmat = [evecf(:,1)*exp(evalf(1)*xcoll(xc)) , evecf(:,2)*exp(evalf(2)*xcoll(xc))];
    Vinvgp{c}(:,xc) = Vmat\(gvalsp{c}(:,xc));
  end
  
  zvalsp{c} = zeros(size(gvalsp{c}));
  for vc = 1:2
    zvalsp{c}(vc,:) = (cumsum(Vinvgp{c}(vc,:)) - 0.5*(Vinvgp{c}(vc,1)+Vinvgp{c}(vc,:)))*dx;
  end
  
  dzvalsp(:,c) = zvalsp{c}(:,end)-zvalsp{c}(:,1);
end

%% solve for boundary conditions , at ends and all junctions
% coefficient order a_1,b_1,a_1',b_1'...a_n,b_n,a_n',a_n',a_(n+1),b_(n+1)
% assume equal phagosome production at all tips and splitting at junctions
% accoding to number of downstream tips

numBranch = length(params.branchpoints);
Mc = zeros(4*numBranch+2);
v = zeros(4*numBranch+2,1);

% branch lens (on main axons)
branchlen = diff([0 params.branchpoints params.L]);
% branch lens (collaterals)
branchlenp = params.branchlen;

% exponential shift to improve scaling
expshift = 0;
shift = 1;

% splitting probabilities to main axon
psplit = params.psplit;

%% equations for joining conditions at branch points
ct = 1; % equation counter
for c = 1:numBranch
  ind = mainind{c};
  
  % first junction equation (for Br)
  % a_i
  Mc(ct,4*(c-1)+1) = evecf(1,1)*exp(evalf(1)*branchlen(c)-expshift);
  % b_i
  Mc(ct,4*(c-1)+2) = evecf(1,2)*exp(evalf(2)*branchlen(c)-expshift);
  % a_i'
  Mc(ct,4*(c-1)+3) = evecf(1,1)*exp(evalf(1)*branchlenp(c)-expshift);
  % b_i'
  Mc(ct,4*(c-1)+4) = evecf(1,2)*exp(evalf(2)*branchlenp(c)-expshift);
  %a_{i+1}
  Mc(ct,4*(c-1)+5) = -evecf(1,1)*shift;
  %b_{i+1}
  Mc(ct,4*(c-1)+6) =  -evecf(1,2)*shift;
  
  % solution vector
  v(ct) =  - evecf(1,:)*dzvals(:,c+1)*shift;
  
  %% junction equation for Ba splitting to main branch
  % a_i
  Mc(ct+1,4*(c-1)+1) = evecf(2,1)*exp(evalf(1)*branchlen(c)-expshift);
  % b_i
  Mc(ct+1,4*(c-1)+2) = evecf(2,2)*exp(evalf(2)*branchlen(c)-expshift);
  % a_{i+1}
  Mc(ct+1,4*(c-1)+5) = -psplit(c)*evecf(2,1)*shift;
  % a_{i+1}
  Mc(ct+1,4*(c-1)+6) = -psplit(c)*evecf(2,2)*shift;
  
  v(ct+1) = -psplit(c)*evecf(2,:)*dzvals(:,c+1)*shift;
  
  %% junction equations for Ba splitting to collateral
  % a_i'
  Mc(ct+2,4*(c-1)+3) = evecf(2,1)*exp(evalf(1)*branchlenp(c)-expshift);
  % b_i'
  Mc(ct+2,4*(c-1)+4) = evecf(2,2)*exp(evalf(2)*branchlenp(c)-expshift);
  % a_{i}
  Mc(ct+2,4*(c-1)+5) =  -(1-psplit(c))*evecf(2,1)*shift;
  % b_{i}
  Mc(ct+2,4*(c-1)+6) = -(1-psplit(c))*evecf(2,2)*shift;
  
  v(ct+2) = -(1-psplit(c))*evecf(2,:)*dzvals(:,c+1)*shift;
  
  %% equation for the tip of the collateral
  % include terms for fusion at the distal tip if necessary
  pfY = params.pft*d.Yt(c);
  %a_i'
  Mc(ct+3,4*(c-1)+3) = (evecf(1,1)-evecf(2,1))*shift;
  % b_i'
  Mc(ct+3,4*(c-1)+4) = (evecf(1,2)-evecf(2,2))*shift;
  v(ct+3) = ((evecf(1,:)-evecf(2,:))*dzvalsp(:,c) + params.kpp(c)*pfY/params.vp + pfY*d.Baup{c}(1))*shift;
  
  ct=ct+4;
end

%% equation for main axon distal tip
pfY = params.pft*d.Yt0;
% a_1
Mc(ct,1) = (evecf(1,1)-evecf(2,1))*shift;
% b_1
Mc(ct,2) = (evecf(1,2)-evecf(2,2))*shift;
v(ct) = ((evecf(1,:)-evecf(2,:))*dzvals(:,1) + params.kpp0*pfY/params.vp + pfY*d.Bau(1))*shift;

% equation for soma
% a_{n+1}
Mc(ct+1,end-1) = evecf(2,1)*exp(evalf(1)*branchlen(end)-expshift);
% b_{n+1}
Mc(ct+1,end) = evecf(2,2)*exp(evalf(2)*branchlen(end)-expshift);
v(ct+1) = 0;

coeff = Mc\v;

%% Ba and Br values on collaterals (and R value at end of each collateral)
Bafp = {}; Brfp = {}; Bsfp = {};
for c = 1:numBranch
  fbs = params.ell*d.Bsup{c}.*(d.kfsap(c).*d.Yap{c} + d.kfsrp(c).*d.Yrp{c});
  
  ind = collind{c};
  xcoll = xvals(ind);
  
  zavals = coeff(4*(c-1)+3:4*(c-1)+4) +zvalsp{c} - zvalsp{c}(:,end);
  
  % bidirectional AVs on collaterals
  Bvecp = evecf*(exp(evalf*xcoll).*zavals);
  Brfp{c} = Bvecp(1,:);
  Bafp{c} = Bvecp(2,:);
  
  Bsfp{c} = (params.khb*(Brfp{c}+Bafp{c}) + fbs)/(params.ks+params.kd + params.kwb);
  Btotp{c} = Brfp{c}+Bafp{c}+Bsfp{c};
  
  fr = params.ell*d.Rup{c}.*(d.kfrap(c).*d.Yap{c} + d.kfrsp(c).*(d.Ysap{c}+d.Ysrp{c}));
  fs = params.kwr/(params.kwr+params.kd)*params.ell*d.Sup{c}.*(d.kfsap(c).*d.Yap{c} + d.kfsrp(c).*d.Yrp{c});
  qvals = (params.ks*(Brfp{c}+Bafp{c}+Bsfp{c}) + fs + fr)/params.vp;
  pvals = (params.kd+params.khr)*ones(size(qvals))/params.vp;
  
  Rfp{c} = solveLinear1ODE(xcoll,pvals,qvals,0);
  
  Rfpend(c) = Rfp{c}(end);
  
  
  Sfp{c} = params.khr*Rfp{c}/(params.kwr+params.kd) +fs/params.kwr;
end
%% Ba and Br values on main axon
zavals(2,:) = 0;
for c = 1:numBranch+1
  ind = mainind{c};
  xshift = xvals(ind) - branchpoints(c);
  
  fbs = params.ell*d.Bsu(ind).*(d.kfsa(c).*d.Ya(ind) + d.kfsr(c).*d.Yr(ind));
  
  zavals = coeff(4*(c-1)+1:4*(c-1)+2) +zvals(:,ind) - zvals(:,ind(end));
  Bvec = evecf*(exp(evalf*xshift).*zavals);
  Brf(ind) = Bvec(1,:);
  Baf(ind) = Bvec(2,:);
  
  Bsf(ind) = (params.khb*(Brf(ind)+Baf(ind)) + fbs)/(params.ks+params.kd+params.kwb);
  
  fr = params.ell*d.Ru(ind).*(d.kfra(c).*d.Ya(ind) + d.kfrs(c).*(d.Ysa(ind)+d.Ysr(ind)));
  fs = params.kwr/(params.kwr+params.kd)*params.ell*d.Su(ind).*(d.kfsa(c).*d.Ya(ind) + d.kfsr(c).*d.Yr(ind));
  qvals = (params.ks*(Brf(ind)+Baf(ind)+Bsf(ind)) + fs + fr)/params.vp;
  pvals = params.kd*ones(size(qvals))/params.vp;
  
  if (c==1)
    Rf(ind) = solveLinear1ODE(xshift,pvals,qvals,0);
  else
    Rf(ind) = solveLinear1ODE(xshift,pvals,qvals,Rfend(c-1)+Rfpend(c-1));
  end
  
  Rfend(c) = Rf(ind(end));
  
  Sf(ind) = params.khr*Rf(ind)/(params.kwr+params.kd) + fs/params.kwr;
end

%% structure with densities of fused but undecayed AVs
fd = struct();
fd.Rf = Rf;
fd.Sf = Sf;
fd.Baf = Baf;
fd.Bsf = Bsf;
fd.Brf = Brf;
if (numBranch>0)
  fd.Rfp = Rfp;
  fd.Sfp = Sfp;
  fd.Bafp = Bafp;
  fd.Brfp = Brfp;
  fd.Bsfp = Bsfp;
else
  fd.Rfp = {};
  fd.Sfp = {};
  fd.Bafp = {};
  fd.Brfp = {};
  fd.Bsfp = {};
end
end