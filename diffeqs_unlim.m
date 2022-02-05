function yp = diffeqs_unlim(x,y,params,Ydata,xdata)

% order of fields:
% Bru, Bau, Ru, Ya, Yr
% order for Ydata: Ya, Yr, Ysa, Ysr

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
    
    Yinterp = interp1(xdata{c},Ydata{c}',x*seglen(c),'linear','extrap');    
    Ya = Yinterp(1); Yr = Yinterp(2); Ysa = Yinterp(3); Ysr = Yinterp(4);
    
    Bru = y(sc+1); Bau = y(sc+2); Ru = y(sc+3);
    
    fs = params.ell*(kfsa*Ya + kfsr*Yr); % fusion with stationary AVs
    Bsu = params.khb*(Bau+Bru)/(params.kwb + params.ks + fs);
    Su = params.khr*Ru/(fs + params.kwr);
    
    
    fr = params.ell*(kfra*Ya + kfrs*(Ysa+Ysr)+kfrr*Yr);
    frb = params.ell*(kfrab*Ya + kfrsb*(Ysa+Ysr)+kfrrb*Yr);
    fab = params.ell*(kfarb*Yr + kfasb*(Ysa+Ysr)+kfaab*Ya);
       
    % dBru/dx
    yp(sc+1) = (-params.khb*Bru + params.pretro*params.kwb*Bsu - (params.ks+frb)*Bru)/params.vpb;
    % dBau/dx
    yp(sc+2) = (params.khb*Bau - (1-params.pretro)*params.kwb*Bsu + (params.ks+fab)*Bau)/params.vpb;
    % dRu/dx
    yp(sc+3) = (params.ks*(Bau + Bru+Bsu) - (params.khr+fr)*Ru + params.kwr*Su)/params.vp;
    
    % rescale equation to segment length
    yp(sc+1:sc+3) = yp(sc+1:sc+3)*seglen(c);
    
    sc=sc+3;
end

end