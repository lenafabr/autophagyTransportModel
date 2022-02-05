function yp = diffeqs_onedone(y,params)

% order of fields:
% Bru, Bau, Ru, Ya, Yr
% first main line segments, then collaterals

% segment lengths
branchlen = diff([0 params.branchpoints params.L]);
branchlenp = params.branchlen;
seglen = [branchlen branchlenp];
nb = length(params.branchlen);

% differential equations on all main segments and collateral branches;
sc=0;
for c = 1:2*nb+1
  
  if (c<=nb+1)
    kfar = params.kfar(c); kfrr = params.kfrr(c); kfra = params.kfra(c); kfsa = params.kfsa(c);
    kfaa = params.kfaa(c); kfas = params.kfas(c); kfrs = params.kfrs(c); kfsr = params.kfsr(c);
    kfarb = params.kfarb(c); kfrrb = params.kfrr(c); kfrab = params.kfra(c);
    kfaab = params.kfaab(c); kfasb = params.kfas(c); kfrsb = params.kfrs(c);
  else
    cc = c-nb-1;
    kfar = params.kfarp(cc); kfrr = params.kfrrp(cc); kfra = params.kfrap(cc); kfsa = params.kfsap(cc);
    kfaa = params.kfaap(cc); kfas = params.kfasp(cc); kfrs = params.kfrsp(cc); kfsr = params.kfsrp(cc);
    kfarb = params.kfarp(cc); kfrrb = params.kfrrp(cc); kfrab = params.kfrap(cc);
    kfaab = params.kfaap(cc); kfasb = params.kfasp(cc); kfrsb = params.kfrsp(cc);
  end
  
  Bru = y(sc+1); Bau = y(sc+2); Ru = y(sc+3); Ya = y(sc+4); Yr = y(sc+5);
  fs = params.ell*(kfsa*Ya + kfsr*Yr);
  Bsu = params.khb*(Bau+Bru)/(params.kwb + params.ks + fs);
  Su = params.khr*Ru/(fs + params.kwr);
  
  
  fyr = params.ell*(kfsr*(Bsu+Su) + kfarb*Bau+kfrrb*Bru+kfrr*Ru);
  fya = params.ell*(kfsa*(Bsu+Su) + kfrab*Bru + kfra*Ru+kfaab*Bau);
  fys = params.ell*(kfas*Bau + kfrsb*Bru + kfrs*Ru);
  
  Ysa = params.khy(1)*Ya/(params.kwy + fys);
  Ysr = params.khy(1)*Yr/(params.kwy + fys);
  
  fr = params.ell*(kfra*Ya + kfrs*(Ysa+Ysr)+kfrr*Yr);
  frb = params.ell*(kfrab*Ya + kfrsb*(Ysa+Ysr)+kfrrb*Yr);
  fab = params.ell*(kfarb*Yr + kfasb*(Ysa+Ysr)+kfaab*Ya);
  
  % dBru/dx
  yp(sc+1) = (-params.khb*Bru + params.pretro*params.kwb*Bsu - (params.ks+frb)*Bru)/params.vpb;
  % dBau/dx
  yp(sc+2) = (params.khb*Bau - (1-params.pretro)*params.kwb*Bsu + (params.ks+fab)*Bau)/params.vpb;
  % dRu/dx
  yp(sc+3) = (params.ks*(Bau + Bru+Bsu) - (params.khr+fr)*Ru + params.kwr*Su)/params.vp;
  % dYa/dx
  yp(sc+4) = ((fya+params.khy(1)).*Ya - params.kwy*Ysa)/params.vya;
  % dYr/dx
  yp(sc+5) = (-(fyr+params.khy(1))*Yr + params.kwy*Ysr)/params.vyr;
  
  % rescale equation to segment length
  yp(sc+1:sc+5) = yp(sc+1:sc+5)*seglen(c);
  
  sc=sc+5;
end