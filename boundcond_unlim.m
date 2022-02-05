function res = boundcond_unlim(ya,yb,params,dens)

% residual of boundary conditions:
% f(ya,yb) = 0 for each boundary condition

% field order: Bru, Bau, Ru, Ya, Yr
% all fields for each main segment first, then each collateral

% conditions at main tip

% tip lysosomes
Yt0 = dens.Yt0;

% condition on Bru(0)
res(1) = (1-params.pft*Yt0)*(params.kpp0 + params.vpb*ya(2)) - params.vpb*ya(1);
% condition on Ru(0)
res(2) = ya(3);
    
% conditions on collateral tips
cc=3;
nb = length(params.branchpoints);
for c= 1:nb
    sc = 3*(nb+c);   
    Yt(c) = dens.Yt(c);
   
    % condition on Brup(0)
    res(cc) = (1-params.pft*Yt(c))*(params.kpp0 + params.vpb*ya(sc+2)) - params.vpb*ya(sc+1);
        
    % condition on Ru(0)
    res(cc+1) = ya(sc+3);    
    
    cc=cc+2;
end

% conditions at each junction
for c = 1:nb
    sc0 = 3*(c-1); %index start for incoming main segment
    sc0p = 3*c; % outgoing main segment
    sc = 3*(nb+c); % collateral
    
    % condition on Bru       
    res(cc) = yb(sc0+1)+yb(sc+1) - ya(sc0p+1);        
    
    % condition on Ru
    res(cc+1) =  yb(sc0+3)+yb(sc+3) - ya(sc0p+3);        
    
    % conditions on Bau    
    res(cc+2) = yb(sc0+2) - params.psplit(c)*ya(sc0p+2);    
    res(cc+3) = yb(sc+2) - (1-params.psplit(c))*ya(sc0p+2);    
    
    cc=cc+4;
end

% condition on Bau(L),soma
sc = 3*nb;
res(cc) = yb(sc+2);


end