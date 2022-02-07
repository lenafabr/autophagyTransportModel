function y = init_guess_unlim(x,params,dens)

% set up initial guess
% field order: Bru, Bau, Ru

sc=0;

nb = length(params.branchpoints);

branchpoints = [0 params.branchpoints params.L];

branchlen = diff([0 params.branchpoints params.L]);
branchlenp = params.branchlen;
seglen = [branchlen branchlenp];

for c = 1:2*nb+1
    
    if (c<=nb+1)
        ind = dens.mainind{c};
        xvals = dens.xvals(ind)-branchpoints(c);     
        y(sc+1) = interp1(xvals,dens.Br(ind),x*seglen(c),'linear','extrap');
        y(sc+2) = interp1(xvals,dens.Ba(ind),x*seglen(c),'linear','extrap');
        y(sc+3) = interp1(xvals,dens.R(ind),x*seglen(c),'linear','extrap');
    else
        cc = c-nb-1;
        xvals = dens.xvals(dens.collind{cc});
        y(sc+1) = interp1(xvals,dens.Brp{cc},x*seglen(c),'linear','extrap');
        y(sc+2) = interp1(xvals,dens.Bap{cc},x*seglen(c),'linear','extrap');
        y(sc+3) = interp1(xvals,dens.Rcoll{cc},x*seglen(c),'linear','extrap');
    end
    
    sc=sc+3;
end

end