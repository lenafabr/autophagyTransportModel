function dens = calcTotalAVdistribution(params,options)
% calculates the total distribution of AVs
% does not include fusion with lysosomes
% WARNING: branchpoints should be offset from xvals

opt = struct();
opt.nx = 1000;
opt.pfsplit = false; % do fusion probabilities depend on radius?

if (exist('options','var'))
  opt = copyStruct(options,opt,1);
end
xvals = linspace(0,params.L,opt.nx);
dx = xvals(2)-xvals(1);

%% first solve for phagosome distribution regardless of fusion

% matrix for solving for anterograde and retrograde bidirectional motion
Bsprefix = params.khb/(params.kwb+params.ks);
% Mrr
M(1,1) = -params.khb + params.pretro*params.kwb*Bsprefix -params.ks;
% Mra
M(1,2) = params.pretro*params.kwb*Bsprefix;
% Mar
M(2,1) = -(1-params.pretro)*params.kwb*Bsprefix;
% Maa
M(2,2) = params.khb - (1-params.pretro)*params.kwb*Bsprefix +params.ks;
M = M/params.vpb;

[evecs,evals] = eig(M);
evals = diag(evals);
%% solve for boundary conditions at ends and all junctions
% coefficient order a_1,b_1,a_1',b_1'...a_n,b_n,a_n',a_n',a_(n+1),b_(n+1)
% assume equal phagosome production at all tips and splitting at junctions
% accoding to number of downstream tips

numBranch = length(params.branchpoints);
Mc = zeros(4*numBranch+2);
v = zeros(4*numBranch+2,1);

% main axon segment lengths
branchlen = diff([0 params.branchpoints params.L]);
% branch lens (collaterals)
branchlenp = params.branchlen;

% exponential shift to improve scaling
expshift = params.L/(numBranch+1)*max(evals);
shift = exp(-expshift);
% splitting probabilities to main axon
psplit = params.psplit;

ct = 1; % equation counter
% equations for joining conditions at branch points
for c = 1:numBranch
  
  % first junction equation (for Br)
  % a_i
  Mc(ct,4*(c-1)+1) = evecs(1,1)*exp(evals(1)*branchlen(c)-expshift);
  % b_i
  Mc(ct,4*(c-1)+2) = evecs(1,2)*exp(evals(2)*branchlen(c)-expshift);
  % a_i'
  Mc(ct,4*(c-1)+3) = evecs(1,1)*exp(evals(1)*branchlenp(c)-expshift);
  % b_i'
  Mc(ct,4*(c-1)+4) = evecs(1,2)*exp(evals(2)*branchlenp(c)-expshift);
  %a_{i+1}
  Mc(ct,4*(c-1)+5) = -evecs(1,1)*shift;
  %b_{i+1}
  Mc(ct,4*(c-1)+6) =  -evecs(1,2)*shift;
  
  % junction equation for Ba splitting to main branch
  % a_i
  Mc(ct+1,4*(c-1)+1) = evecs(2,1)*exp(evals(1)*branchlen(c)-expshift);
  % b_i
  Mc(ct+1,4*(c-1)+2) = evecs(2,2)*exp(evals(2)*branchlen(c)-expshift);
  % a_{i+1}
  Mc(ct+1,4*(c-1)+5) = -psplit(c)*evecs(2,1)*shift;
  % a_{i+1}
  Mc(ct+1,4*(c-1)+6) = -psplit(c)*evecs(2,2)*shift;
  
  % junction equations for Ba splitting to collateral
  % a_i'
  Mc(ct+2,4*(c-1)+3) = evecs(2,1)*exp(evals(1)*branchlenp(c)-expshift);
  % b_i'
  Mc(ct+2,4*(c-1)+4) = evecs(2,2)*exp(evals(2)*branchlenp(c)-expshift);
  % a_{i}
  Mc(ct+2,4*(c-1)+5) =  -(1-psplit(c))*evecs(2,1)*shift;
  % b_{i}
  Mc(ct+2,4*(c-1)+6) = -(1-psplit(c))*evecs(2,2)*shift;
  
  % equation for the tip of the collateral
  %a_i'
  Mc(ct+3,4*(c-1)+3) = (evecs(1,1)-evecs(2,1))*shift;
  % b_i'
  Mc(ct+3,4*(c-1)+4) = (evecs(1,2)-evecs(2,2))*shift;
  
  % solution vector
  v(ct+3) = params.kpp(c)/params.vpb*shift;
  
  ct=ct+4;
end

% equation for main axon tip
% a_1
Mc(ct,1) = (evecs(1,1) - evecs(2,1))*shift;
% b_1
Mc(ct,2) = (evecs(1,2)-evecs(2,2))*shift;
v(ct) = params.kpp0/params.vpb*shift;

% equation for soma
% a_{n+1}
Mc(ct+1,end-1) = evecs(2,1)*exp(evals(1)*branchlen(end)-expshift);
% b_{n+1}
Mc(ct+1,end) = evecs(2,2)*exp(evals(2)*branchlen(end)-expshift);

coeff = Mc\v;

branchpoints = [0 params.branchpoints params.L];

% Ba and Br values on collaterals (and R value at end of each collateral)
Bap = cell(1,numBranch); Brp = cell(1,numBranch); Bsp = cell(1,numBranch); Sp = cell(1,numBranch);
collind = cell(1,numBranch); Rcoll = cell(1,numBranch); Btotp = cell(1,numBranch);

Rcollend = zeros(1,numBranch);

for c = 1:numBranch
  ind = find(xvals<=params.branchlen(c));
  collind{c} = ind;
  
  % bidirectional AV on collaterals
  Bvecp = coeff(4*(c-1)+3)*evecs(:,1)*exp(evals(1)*xvals(ind)) + coeff(4*(c-1)+4)*evecs(:,2)*exp(evals(2)*xvals(ind));
  Brp{c} = Bvecp(1,:);
  Bap{c} = Bvecp(2,:);
  Bsp{c} = Bsprefix*(Brp{c}+Bap{c});
  Btotp{c} = Brp{c}+Bap{c}+Bsp{c};
  
  
  % Ru on collateral
  Rcoll{c} = params.ks/params.vp*(cumsum(Btotp{c})-0.5*(Btotp{c}(1)+Btotp{c}))*dx;
  %totsum = (sum(Btotp{c}) - 0.5*(Btotp{c}(1)+Btotp{c}(end)))*dx;
  Rcollend(c) = Rcoll{c}(end);
  
  Sp{c} = params.khr/params.kwr*Rcoll{c};
  %params.ks/params.vp*totsum;
end

% Ba and Br values on main axon
Br = zeros(size(xvals));
Ba = Br;
Bs = Br;
R = Br;
S = Br;
Btot = Br;

Rend = zeros(1,numBranch+1);
mainind = cell(1,numBranch+1);
for c = 1:numBranch+1
  ind = find(xvals>=branchpoints(c) & xvals<=branchpoints(c+1));
  mainind{c} = ind;
  
  xshift = xvals(ind) - branchpoints(c);
  
  % on main axon
  Bvec = coeff(4*(c-1)+1)*evecs(:,1)*exp(evals(1)*xshift) + coeff(4*(c-1)+2)*evecs(:,2)*exp(evals(2)*xshift);
  Br(ind) = Bvec(1,:);
  Ba(ind) = Bvec(2,:);
  Bs(ind) = Bsprefix*(Br(ind)+Ba(ind));
  Btot(ind) = Br(ind)+Ba(ind)+Bs(ind);
  
  % integrate over this main branch only
  R(ind) = (cumsum(Btot(ind))-0.5*(Btot(ind(1))+Btot(ind)))*(dx*params.ks/params.vp);
  
  if (c>1)
    R(ind) = R(ind)+Rend(c-1)+Rcollend(c-1);
  end
  
  Rend(c) = R(ind(end)); % retrograde phagos at end of this branch
  
  % fully retrograde ones that went stationary
  S(ind) = params.khr/params.kwr*R(ind);
end

%% output structure
dens = struct();
dens.xvals = xvals;
dens.R = R; % fully retrograde AVs on main branch
dens.S = S; % fully retrograde AVs that went stationary
dens.Br = Br; % bidirectional retrograde, main branch
dens.Ba = Ba; % bidirectiona anterograde, main branch
dens.Bs = Bs; % bidirectional stationary main branch

% corresponding quantities on collaterals
dens.Rcoll = Rcoll;
dens.Sp = Sp;
dens.Brp = Brp;
dens.Bap = Bap;
dens.Bsp = Bsp;

dens.mainind = mainind;
dens.collind = collind;

end