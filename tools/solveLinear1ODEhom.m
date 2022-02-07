function yvals = solveLinear1ODEhom(xvals,pvals,yfix,fixy0)
% solve homogeneous 1st order ODE y' = p(x)*y
% if fixy0=true then y(0) = yfix
% if fixy0 = false then y(L) = yfix

dx = xvals(2)-xvals(1);
pint = (cumsum(pvals)-0.5*(pvals(1)+pvals))*dx;

if (fixy0)
    yvals = yfix*exp(pint);
else
    yvals = yfix*exp(-pint(end) +pint);
end

end