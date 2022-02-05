function res = boundcond_onedone(ya,yb,params)

% residual of boundary conditions:
% f(ya,yb) = 0 for each boundary condition

% field order: Bru, Bau, Ru, Ya, Yr
% all fields for each main segment first, then each collateral

% conditions at main tip

% tip lysosomes
Yt0 = params.vya*ya(4)/(params.kye + params.pft*(params.vpb*ya(2) + params.kpp0));
% condition on Bru(0)
res(1) = (1-params.pft*Yt0)*(params.kpp0 + params.vpb*ya(2)) - params.vpb*ya(1);
% condition on Ru(0)
res(2) = ya(3);
% condition on Yr(0)
res(3) = params.kye*Yt0 - params.vyr*ya(5);
    
% conditions on collateral tips
cc=4;
nb = length(params.branchpoints);
for c= 1:nb
    sc = 5*(nb+c);   
    Yt(c) = params.vya*ya(sc+4)/(params.kye + params.pft*(params.vpb*ya(sc+2) + params.kpp(c)));
    
    % condition on Brup(0)
    res(cc) = (1-params.pft*Yt(c))*(params.kpp0 + params.vpb*ya(sc+2)) - params.vpb*ya(sc+1);
        
    % condition on Ru(0)
    res(cc+1) = ya(sc+3);
    
    % condition on Yr(0)    
    res(cc+2) = params.kye*Yt(c) - params.vyr*ya(sc+5);
    
    cc=cc+3;
end

% conditions at each junction
for c = 1:nb
    sc0 = 5*(c-1); %index start for incoming main segment
    sc0p = 5*c; % outgoing main segment
    sc = 5*(nb+c); % collateral
    
    % condition on Bru       
    res(cc) = yb(sc0+1)+yb(sc+1) - ya(sc0p+1);        
    
    % condition on Ru
    res(cc+1) =  yb(sc0+3)+yb(sc+3) - ya(sc0p+3);        
    
    % condition on Yr    
    res(cc+2) = yb(sc0+5)+yb(sc+5) - ya(sc0p+5);
    
    % conditions on Ba    
    res(cc+3) = yb(sc0+2) - params.psplit(c)*ya(sc0p+2);    
    res(cc+4) = yb(sc+2) - (1-params.psplit(c))*ya(sc0p+2);
    
    % conditions on Ya
    res(cc+5) = yb(sc0+4) - params.psplit(c)*ya(sc0p+4); 
    res(cc+6) = yb(sc+4) - (1-params.psplit(c))*ya(sc0p+4);
    
    cc=cc+7;
end

% condition on Bau(L),soma
sc = 5*nb;
res(cc) = yb(sc+2);

% condition on Ya(L), soma
res(cc+1) = params.kpy - params.vya*yb(sc+4);

end