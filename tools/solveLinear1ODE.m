function yvals = solveLinear1ODE(xvals,pvals,qvals,y0)
% solve a linear first-order differential equation using an integration
% factor with numerical integration
% y' + p(x) y = q(x)
% y(0) = y0
% assumes xvals are evenly spaced

dx = xvals(2)-xvals(1);
% integral of P 
intP = (cumsum(pvals)-0.5*(pvals(1)+pvals))*dx;
for xc = 1:length(xvals)
    integ = exp(intP(1:xc)-intP(xc)).*qvals(1:xc);       
   % assume last point is equal to the previous one
    yvals(xc) = (sum(integ)-0.5*integ(1)-0.5*integ(xc))*dx + y0*exp(-intP(xc));
end
