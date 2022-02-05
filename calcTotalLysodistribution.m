function dens = calcTotalLysodistribution(params,dens,options)
% solves for distributions of lysosomes
% total AV densities precalculated with calcTotalAVdistribution and stored in dens

opt = struct();
opt.dotipfusion = true;
if (exist('options','var'))
  opt = copyStruct(options,opt,'addnew',1);
end

Yap = {}; Yrp = {}; Ysap = {}; Ysrp = {};

branchpoints = [0 params.branchpoints params.L];
numBranch = length(params.branchpoints);
xvals = dens.xvals;
psplit = params.psplit;

%% solve for lysosomes
Yt0 = 0; Yt = zeros(numBranch,1);

for c=numBranch+1 : -1 : 1
  
  % anterograde lysosomes on main section
  ind= dens.mainind{c};
  xshift = xvals(ind) - branchpoints(c);
  
  % fusing with antero lyso
  fuserate = params.ell*(params.kfra(c)*dens.R(ind)+params.kfrab(c)*dens.Br(ind) + params.kfaab(c)*dens.Ba(ind)+ params.kfsa(c)*(dens.S(ind)+dens.Bs(ind)));
  
  % fusing with stationary lyso
  fuserateS = params.ell*(params.kfrs(c)*dens.R(ind)+params.kfrsb(c)*dens.Br(ind) + params.kfasb(c)*dens.Ba(ind));
  
  qvals = (fuserate + params.kda + params.khy.*(1 - params.kwy./(fuserateS + params.kwy+params.kds)))/params.vya;
  
  if (c==numBranch+1) % most proximal region (adjacent to soma)
    Ya(ind) = solveLinear1ODEhom(xshift,qvals,params.kpy/params.vya,false);
  else % more downstream regions
    Ya(ind) = solveLinear1ODEhom(xshift,qvals,Yajunc(c+1)*psplit(c),false);
  end
  Yajunc(c) = Ya(ind(1)); % value at the junction on the distal side of this segment
  
  Ysa(ind) = params.khy.*Ya(ind)./(fuserateS + params.kwy+params.kds);
  
  %  lysosomes on collateral branch
  if (c<numBranch+1)
    ind = dens.collind{c};
    xcoll = xvals(ind);
    
    % anterograde
    fuserate = params.ell*(params.kfrap(c)*dens.Rcoll{c}+params.kfrapb(c)*dens.Brp{c} + params.kfaapb(c)*dens.Bap{c}+params.kfsap(c)*(dens.Sp{c}+dens.Bsp{c}));
    fuserateS = params.ell*(params.kfrsp(c)*dens.Rcoll{c}+params.kfrspb(c)*dens.Brp{c} + params.kfaspb(c)*dens.Bap{c});
    
    qvals = (fuserate + params.kda + params.khy.*(1 - params.kwy./(fuserateS + params.kwy)))/params.vya;
    Yap{c} = solveLinear1ODEhom(xcoll,qvals,Yajunc(c+1)*(1-psplit(c)),false);
    Ysap{c} = params.khy.*Yap{c}./(fuserateS + params.kwy+params.kds);
    
    if (opt.dotipfusion)
      % tip lysosomes on collateral
      
      if (isinf(params.kpe)) % phagosomes bounce off distal end
        Yt(c) = params.vya*Yap{c}(1)/(params.kye + params.pft*(params.vp*dens.Bap{c}(1)+params.kpp(c)));
      else % phagosomes pile up at distal end
        error('not set up')
        %Yt(c) = params.vya*Yap{c}(1)/(params.kye + params.kft*Pt(c));
      end
      
    end
    
    % retrograde on collateral
    
    % fusing with retro lyso
    fuserate = params.ell*(params.kfrrp(c)*dens.Rcoll{c}+params.kfrrpb(c)*dens.Brp{c}+params.kfarpb(c)*dens.Bap{c}+params.kfsrp(c)*(dens.Sp{c}+dens.Bsp{c}));
    qvals = (-fuserate -params.kdr + params.khy.*(-1 + params.kwy./(fuserateS + params.kwy)))/params.vyr;
    
    if (opt.dotipfusion)
      Yrbound = params.kye*Yt(c)/params.vyr;
    else
      Yrbound = Yap{c}(1)*params.vya/params.vyr;
    end
    
    Yrp{c} = solveLinear1ODEhom(xcoll,qvals,Yrbound,true);
    
    Ysrp{c} = params.khy.*Yrp{c}./(fuserateS + params.kwy+params.kds);
  end
end

% tip lysosomes on main branch
if (opt.dotipfusion)
  % main axon tip
  if (isinf(params.kpe)) % phagosomes bounce off distal end
    Yt0 = params.vya*Ya(1)/(params.kye + params.pft*(params.vp*dens.Ba(1)+params.kpp0));
  else % phagosomes pile up at distal end
    Yt0 = params.vya*Ya(1)/(params.kye + params.kft*Pt0);
  end
end

Yr = zeros(1,opt.nx); Ysr = Yr;
Yrend = zeros(1,numBranch+1);

% go upstream to solve for retrograde lysos on main axon
for c = 1:numBranch+1
  ind = dens.mainind{c};
  xshift = xvals(ind) - branchpoints(c);
  
  fuserate = params.ell*(params.kfrr(c)*dens.R(ind)+params.kfrrb(c)*dens.Br(ind) + params.kfarb(c)*dens.Ba(ind) + params.kfsr(c)*(dens.S(ind)+dens.Bs(ind)));
  fuserateS = params.ell*(params.kfrs(c)*dens.R(ind)+params.kfrsb(c)*dens.Br(ind) + params.kfasb(c)*dens.Ba(ind));
  
  qvals = (-fuserate -params.kdr + params.khy.*(-1 + params.kwy./(fuserateS + params.kwy)))/params.vyr;
  
  if (c==1)
    if (opt.dotipfusion)
      Yrbound = params.kye*Yt0/params.vyr;
    else
      Yrbound = Ya(1)*params.vya/params.vyr;
    end
    
    Yr(ind) = solveLinear1ODEhom(xshift,qvals,Yrbound,true);
  else
    Yr(ind) = solveLinear1ODEhom(xshift,qvals,Yrp{c-1}(end)+Yrend(c-1),true);
  end
  Yrend(c) = Yr(ind(end));
  
  Ysr(ind) = params.khy.*Yr(ind)./(fuserateS + params.kwy+params.kds);
  
end

%% add to density structure
dens.Ya = Ya;
dens.Yr = Yr;
dens.Yap = Yap;
dens.Yrp = Yrp;
dens.Ysr = Ysr;
dens.Ysa = Ysa;
dens.Ysrp = Ysrp;
dens.Ysap = Ysap;
dens.Yt0 = Yt0;
dens.Yt = Yt;

end