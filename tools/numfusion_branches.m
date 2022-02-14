function fd = numfusion_branches(params,d)
% given output densities from bidirectionAV_branches
% calculate the density of fusions (ie: # fusion * AV density) in main axon
% and branches
% d = intput densities structure
% returns structure containing all densities
% interpolate along segments to match up boundary conditions better


xvals = d.xvals;
dx = xvals(2)-xvals(1);
nb = length(params.branchpoints);
branchpoints = [0 params.branchpoints params.L];
% branch lens (on main axons)
branchlen = diff([0 params.branchpoints params.L]);
% branch lens (collaterals)
branchlenp = params.branchlen;

%% on main axon

Mf(1,1) = -(params.ks + params.khb)+ params.khb*params.kwb*params.pretro/(params.ks + params.kwb);
Mf(1,2) = params.khb*params.kwb*params.pretro/(params.ks + params.kwb);
Mf(2,1) = - params.khb*params.kwb*(1-params.pretro)/(params.ks + params.kwb);
Mf(2,2) = (params.ks + params.khb)- params.khb*params.kwb*(1-params.pretro)/(params.ks + params.kwb);
Mf = Mf/params.vp;

[evecf,evalf] = eig(Mf);
evalf = diag(evalf);

% sort eigenvalues (2nd one larger);
[a,b] = sort(evalf);
evecf = evecf(:,b);
%% solve for z and w functions on main axon
mainind = {};
for c = 1:nb+1
    ind = find(xvals>=branchpoints(c) & xvals<=branchpoints(c+1));
    mainind{c} = ind;
    xshift = xvals(ind) - branchpoints(c);
    
    fbr = params.ell*d.Br(ind).*(d.kfra(c).*d.Ya(ind) + d.kfrs(c).*(d.Ysa(ind)+d.Ysr(ind)) + d.kfrr(c).*d.Yr(ind));
    fbs = params.ell*d.Bs(ind).*(d.kfsa(c).*d.Ya(ind) + d.kfsr(c).*d.Yr(ind));
    fba = params.ell*d.Ba(ind).*(d.kfar(c).*d.Yr(ind) + d.kfas(c).*(d.Ysa(ind)+d.Ysr(ind))+ d.kfaa(c).*d.Ya(ind));
    gvals = [fbr + params.kwb*params.pretro/(params.ks + params.kwb)*fbs; ...
            -fba - params.kwb*(1-params.pretro)/(params.ks+params.kwb)*fbs];
    gvals = gvals/params.vp;

    
    Vinvg = zeros(2,length(ind));
    for xc = 1:length(xshift)
        Vmat = [evecf(:,1)*exp(evalf(1)*xshift(xc)) , evecf(:,2)];
        % solve for z-prime-tilde
        Vinvg(:,xc) = Vmat\gvals(:,xc);           
    end   
      
    %z_1
    zwvals(1,ind) = (cumsum(Vinvg(1,:)) - 0.5*(Vinvg(1,1)+Vinvg(1,:)))*dx;    
    %w_2    
    uvals = fliplr(branchlen(c) - xshift);
    pvals = -evalf(2)*ones(size(xshift));
    qvals = fliplr(Vinvg(2,:));
    zwvals(2,ind) = solveLinear1ODE(uvals,-pvals,-qvals,0);
    zwvals(2,ind) = fliplr(zwvals(2,ind));
    
    % save z_1(L)*exp(lam1*L)
    zendexp(c) = zwvals(1,ind(end))*exp(evalf(1)*branchlen(c));
    % w2(0)
    wstart(c) = zwvals(2,ind(1));
end  


%% along collateral branches
% v1*exp(l1*x)*z1 + v2*exp(l2*x)*z2 evaluated at end of collateral

zvalsp = {}; zendexpP = [];
for c = 1:nb
    %%
    collind{c} = find(xvals<=params.branchlen(c));
    ind = collind{c};
    xcoll = xvals(ind);
    
    fbr = params.ell*d.Brp{c}.*(d.kfrap(c).*d.Yap{c} + d.kfrsp(c).*(d.Ysap{c}+d.Ysrp{c}) + d.kfrrp(c).*d.Yrp{c});
    fbs = params.ell*d.Bsp{c}.*(d.kfsap(c).*d.Yap{c} + d.kfsrp(c).*d.Yrp{c});
    fba = params.ell*d.Bap{c}.*(d.kfarp(c).*d.Yrp{c} + d.kfasp(c).*(d.Ysap{c}+d.Ysrp{c})+d.kfaap(c).*d.Yap{c});
    
    %%
    gvalsp{c} = [fbr + params.kwb*params.pretro/(params.ks + params.kwb)*fbs; ...
        -fba - params.kwb*(1-params.pretro)/(params.ks+params.kwb)*fbs];
    gvalsp{c} = gvalsp{c}/params.vp;
    
    Vinvgp{c} = zeros(size(gvalsp{c}));
    for xc = 1:length(xcoll)
        Vmat = [evecf(:,1)*exp(evalf(1)*xcoll(xc)) , evecf(:,2)];
        Vinvgp{c}(:,xc) = Vmat\(gvalsp{c}(:,xc));
    end       
    
    zwvalsp{c} = zeros(size(gvalsp{c}));
    %z_1
    zwvalsp{c}(1,:) = (cumsum(Vinvgp{c}(1,:)) - 0.5*(Vinvgp{c}(1,1)+Vinvgp{c}(1,:)))*dx;
    % w_2
    uvals = fliplr(branchlenp(c)-xcoll);
    pvals = -evalf(2)*ones(size(xcoll));
    qvals = fliplr(Vinvgp{c}(2,:));
    zwvalsp{c}(2,:) = solveLinear1ODE(uvals,-pvals,-qvals,0);
    zwvalsp{c}(2,:) = fliplr(zwvalsp{c}(2,:));
    
    % save z_1(L)*exp(lam1*L)
    zendexpP(c) = zwvalsp{c}(1,end)*exp(evalf(1)*branchlenp(c));
    wstartP(c) = zwvalsp{c}(2,1); 
end

%% solve for boundary conditions for number of fusions, at ends and all junctions
% coefficient order a_1,b_1,a_1',b_1'...a_n,b_n,a_n',a_n',a_(n+1),b_(n+1)
% assume equal phagosome production at all tips and splitting at junctions
% accoding to number of downstream tips

nb = length(params.branchpoints);
Mc = zeros(4*nb+2);
v = zeros(4*nb+2,1);

% branch lens (on main axons)
branchlen = diff([0 params.branchpoints params.L]);
% branch lens (collaterals)
branchlenp = params.branchlen;


% splitting probabilities to main axon
psplit = params.psplit;

%psplit = 0.5*ones(1,nb);

% production rates at each junction tip
kprod = params.kpp;
% production at tip of main axon;
kprod0 = params.kpp;

%% equations for joining conditions at branch points
ct = 1; % equation counter
for c = 1:nb
    ind = mainind{c};        
    
    % first junction equation (for Br)
    % a1_i
    Mc(ct,4*(c-1)+1) = evecf(1,1)*exp(evalf(1)*branchlen(c));
    % b2_i
    Mc(ct,4*(c-1)+2) = evecf(1,2);
    % a1_i'
    Mc(ct,4*(c-1)+3) = evecf(1,1)*exp(evalf(1)*branchlenp(c));
    % b2_i'
    Mc(ct,4*(c-1)+4) = evecf(1,2);
    %a1_{i+1}
    Mc(ct,4*(c-1)+5) = -evecf(1,1);
    %b2_{i+1}
    Mc(ct,4*(c-1)+6) =  -evecf(1,2)*exp(-evalf(2)*branchlen(c+1));
    
    % solution vector
    v(ct) =  - evecf(1,1)*(zendexp(c) + zendexpP(c)) + evecf(1,2)*wstart(c+1);    
    
    %% junction equation for Ba splitting to main branch
    % a1_i
    Mc(ct+1,4*(c-1)+1) = evecf(2,1)*exp(evalf(1)*branchlen(c));
    % b2_i
    Mc(ct+1,4*(c-1)+2) = evecf(2,2);
    % a1_{i+1}
    Mc(ct+1,4*(c-1)+5) = -psplit(c)*evecf(2,1);
    % b2_{i+1}
    Mc(ct+1,4*(c-1)+6) = -psplit(c)*evecf(2,2)*exp(-evalf(2)*branchlen(c+1));
    
    v(ct+1) = psplit(c)*evecf(2,2)*wstart(c+1) - zendexp(c)*evecf(2,1);
    
    %% junction equations for Ba splitting to collateral
    % a1_i'
    Mc(ct+2,4*(c-1)+3) = evecf(2,1)*exp(evalf(1)*branchlenp(c));    
    % b2_i'
    Mc(ct+2,4*(c-1)+4) = evecf(2,2);    
    % a1_{i+1}
    Mc(ct+2,4*(c-1)+5) =  -(1-psplit(c))*evecf(2,1);
    % b2_{i+1}
    Mc(ct+2,4*(c-1)+6) = -(1-psplit(c))*evecf(2,2)*exp(-evalf(2)*branchlen(c+1));
    
    v(ct+2) = (1-psplit(c))*evecf(2,2)*wstart(c+1) - zendexpP(c)*evecf(2,1);
    
    %% equation for the tip of the collateral
    % include terms for fusion at the distal tip if necessary
    pfY = params.pft*d.Yt(c);
    %a1_i'
    Mc(ct+3,4*(c-1)+3) = (evecf(1,1)-evecf(2,1));    
    % b2_i'
    Mc(ct+3,4*(c-1)+4) = (evecf(1,2)-evecf(2,2))*exp(-evalf(2)*branchlenp(c));    
    v(ct+3) = (evecf(2,2)-evecf(1,2))*wstartP(c) + (params.kpp(c)/params.vp+d.Bap{c}(1))*pfY;    
    
    ct=ct+4;
end
    
%% equation for main axon distal tip
pfY = params.pft*d.Yt0;
% a1_1
Mc(ct,1) = (evecf(1,1)-evecf(2,1));
% b2_1
Mc(ct,2) = (evecf(1,2)-evecf(2,2))*exp(-evalf(2)*branchlen(1));
v(ct) = (evecf(2,2)-evecf(1,2))*wstart(1) + (params.kpp0/params.vp+d.Ba(1))*pfY;

% equation for soma
% a1_{n+1}
Mc(ct+1,end-1) = evecf(2,1)*exp(evalf(1)*branchlen(end));
% b2_{n+1}
Mc(ct+1,end) = evecf(2,2);
v(ct+1) = -zendexp(end)*evecf(2,1);

coeff = Mc\v;

%% Ba and Br values on collaterals (and R value at end of each collateral)
Bafp = {}; Brfp = {}; Bsfp = {}; 
for c = 1:nb
    fbs = params.ell*d.Bsp{c}.*(d.kfsap(c).*d.Yap{c} + d.kfsrp(c).*d.Yrp{c});
    
    ind = collind{c};
    xcoll = xvals(ind);
    
    % bidirectional AVs on collaterals
    expvals = [exp(evalf(1)*xcoll); exp(evalf(2)*(xcoll-branchlenp(c)))];    
    inhomvals = [zwvalsp{c}(1,:).*expvals(1,:); zwvalsp{c}(2,:)];
    Bvecp = evecf*(expvals.*coeff(4*(c-1)+3:4*(c-1)+4) + inhomvals);
    
    %Bvec = evecf*(exp(evalf*xshift).*zavals);        
    % z(x) + a = z(x) - z(L) + alpha
    %zavals = coeff(4*(c-1)+3:4*(c-1)+4) +zvalsp{c} - zvalsp{c}(:,end);
    
   
    %Bvecp = evecf*(exp(evalf*xcoll).*zavals);    
    Brfp{c} = Bvecp(1,:);
    Bafp{c} = Bvecp(2,:); 
  
    Bsfp{c} = (params.khb*(Brfp{c}+Bafp{c}) + fbs)/(params.ks+params.kwb);    
    Btotp{c} = Brfp{c}+Bafp{c}+Bsfp{c};
    
    fr = params.ell*d.Rcoll{c}.*(d.kfrap(c).*d.Yap{c} + d.kfrsp(c).*(d.Ysap{c}+d.Ysrp{c})+d.kfrrp(c).*d.Yrp{c});
    fs = params.ell*d.Sp{c}.*(d.kfsap(c).*d.Yap{c} + d.kfsrp(c).*d.Yrp{c});    
    qvals = (params.ks*(Brfp{c}+Bafp{c}+Bsfp{c}) + fs + fr)/params.vp;
    
    Rfp{c} = (cumsum(qvals)-0.5*(qvals(1)+qvals))*dx;            
    Rfpend(c) = Rfp{c}(end);
    
    Sfp{c} = (params.khr*Rfp{c}+fs)/params.kwr;
end
    
    
%     
%     % Ru on collateral
%     Rcoll{c} = params.ks/params.vp*(cumsum(Btotp{c})-0.5*(Btotp{c}(1)+Btotp{c}))*dx;    
%     %totsum = (sum(Btotp{c}) - 0.5*(Btotp{c}(1)+Btotp{c}(end)))*dx;
%     Rcollend(c) = Rcoll{c}(end);
%     
%     Sp{c} = params.khr/params.kwr*Rcoll{c};
%     %params.ks/params.vp*totsum;

%% Ba and Br values on main axon
for c = 1:nb+1
    ind = mainind{c};    
    xshift = xvals(ind) - branchpoints(c);
    
    fbs = params.ell*d.Bs(ind).*(d.kfsa(c).*d.Ya(ind) + d.kfsr(c).*d.Yr(ind));
        
    expvals = [exp(evalf(1)*xshift); exp(evalf(2)*(xshift-branchlen(c)))];    
    inhomvals = [zwvals(1,ind).*expvals(1,:); zwvals(2,ind)];
    Bvec = evecf*(expvals.*coeff(4*(c-1)+1:4*(c-1)+2) + inhomvals);
    %Bvec = evecf*(exp(evalf*xshift).*zavals);
    Brf(ind) = Bvec(1,:);
    Baf(ind) = Bvec(2,:); 
  
    Bsf(ind) = (params.khb*(Brf(ind)+Baf(ind)) + fbs)/(params.ks+params.kwb);
    %Bsfp{c} = (params.khb*(Brfp{c}+Bafp{c}) + fbs)/(params.ks+params.kwb);
    %Btotp{c} = Brfp{c}+Bafp{c}+Bsfp{c};
    
    fr = params.ell*d.R(ind).*(d.kfra(c).*d.Ya(ind) + d.kfrs(c).*(d.Ysa(ind)+d.Ysr(ind))+d.kfrr(c).*d.Yr(ind));
    fs = params.ell*d.S(ind).*(d.kfsa(c).*d.Ya(ind) + d.kfsr(c).*d.Yr(ind));    
    qvals = (params.ks*(Brf(ind)+Baf(ind)+Bsf(ind)) + fs + fr)/params.vp;
    
    Rf(ind) = (cumsum(qvals)-0.5*(qvals(1)+qvals))*dx;
        
    if (c>1)
        Rf(ind) = Rf(ind)+Rfend(c-1)+Rfpend(c-1);
    end
    
    Rfend(c) = Rf(ind(end));
    
    Sf(ind) = (params.khr*Rf(ind)+fs)/params.kwr;
end

% %% check tip boundary conditions
% c=1;
% %Brfp{c}(1)
% Bafp{c}(1) + (params.kpp(c)/params.vp + d.Bap{c}(1))*params.pft*d.Yt(c)
% 
% %% check derivatives
% c=1;
% fba = params.ell*d.Bap{c}.*(d.kfarp(c).*d.Yrp{c} + d.kfasp(c).*(d.Ysap{c}+d.Ysrp{c})+d.kfaap(c).*d.Yap{c});
% xcoll = xvals(collind{1});
% dx = xcoll(2)-xcoll(1);
% dBa = (Bafp{c}(3:end)-Bafp{c}(1:end-2))/(2*dx);
% dBarhs = (params.ks + params.khb)*Bafp{c} - params.kwb*(1-params.pretro)*Bsfp{c} - fba;
% plot(xcoll(2:end-1),dBa,xcoll(2:end-1), dBarhs(2:end-1))
% 
%%
% dBr = (Brf(3:end)-Brf(1:end-2))/(2*dx);
% dBrrhs = -(params.ks + params.khb)*Brf + params.kwb*params.pretro*Bsf + fbr;
% plot(xshift(2:end-1),dBr,xshift(2:end-1), dBrrhs(2:end-1))
%% structure with densities of fusions
fd = struct();
fd.Rf = Rf;
fd.Sf = Sf;
fd.Baf = Baf;
fd.Bsf = Bsf;
fd.Brf = Brf;
fd.mainind = mainind;
if (nb>0)
fd.Rfp = Rfp;
fd.Sfp = Sfp;
fd.Bafp = Bafp;
fd.Brfp = Brfp;
fd.Bsfp = Bsfp;
fd.collind = collind;
else
    fd.Rfp = {};
    fd.Sfp = {};
    fd.Bafp = {};
    fd.Brfp = {};
    fd.Bsfp = {};
    fd.collind = {};
end
end