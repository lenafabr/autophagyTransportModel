function y = init_guess_onedone(params)

% set up initial guess
% field order: Bru, Bau, Ru, Ya, Yr

sc=0;
for c = 1:2*length(params.branchpoints)+1
  y(sc+1) = params.kpp0/params.vp;
  y(sc+2) = y(1);
  y(sc+3) = params.kpp0/params.vp;
  y(sc+4) = params.kpy/params.vya;
  y(sc+5) = params.kpy/params.vyr;
  
  sc=sc+5;
end

end