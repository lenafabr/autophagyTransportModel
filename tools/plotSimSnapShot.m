function plotSimSnapShot(snaptime,grouplist,tvals)
%% PLOTSIMSNAPSHOT prints a snapshot from agent-based simulations
%% at a given time snaptime
%% grouplist is a partGroupObj generated using readsnapshot.m
  [lyso_img,~,lyso_alpha] = imread('lyso_img.png'); lyso_img = flipud(lyso_img); lyso_alpha = flipud(lyso_alpha);
  [SU_img,~,SU_alpha] = imread('SU_img.png'); SU_img = flipud(SU_img); SU_alpha = flipud(SU_alpha);
  [SF_img,~,SF_alpha] = imread('SF_img.png'); SF_img = flipud(SF_img); SF_alpha = flipud(SF_alpha);
  [RU_img,~,RU_alpha] = imread('RU_img.png'); RU_img = flipud(RU_img); RU_alpha = flipud(RU_alpha);
  [RF_img,~,RF_alpha] = imread('RF_img.png'); RF_img = flipud(RF_img); RF_alpha = flipud(RF_alpha);
  
  snapind = find(tvals>snaptime,1);
  maxnpart = max(grouplist.npart);
  yinds = zeros(maxnpart,1);

  Lreal = 1055;
  vreal = 0.75;
  tau = Lreal/vreal;

  sizescl = 10;
  pbscl = 6;
  xmax = Lreal;
  ymax = maxnpart;
  xscl = pbscl*ymax/xmax;

  clf
  cursnapshot = grouplist.snapshot(snapind);
  pidvals = cursnapshot.id;

  delinds = ~ismember(yinds,pidvals);
  yinds(delinds) = 0;

  newinds = ~ismember(pidvals,yinds); %pid indices that are not in yinds
  if(any(newinds))
    reminds = find(yinds==0);
    if(nnz(reminds)==1)
      assigninds = reminds;
    else
      assigninds = randsample(reminds,nnz(newinds));
    end
    yinds(assigninds) = pidvals(newinds);
  end
  [~,ia,ib] = intersect(pidvals,yinds,'stable');

  lysoinds = find(cursnapshot.type==1);
  SUinds = find(cursnapshot.type==2 & cursnapshot.state==1 & cursnapshot.nfuse==0);
  SFinds = find(cursnapshot.type==2 & cursnapshot.state==1 & cursnapshot.nfuse>0);
  RUinds = find(cursnapshot.type==2 & cursnapshot.state==2 & cursnapshot.nfuse==0);
  RFinds = find(cursnapshot.type==2 & cursnapshot.state==2 & cursnapshot.nfuse>0);

  allpos = Lreal*(1-cursnapshot.pos);

  hold on
  set(gca,'Color','white')
  if(~isempty(lysoinds))
    xvals = allpos(lysoinds)+[0,sizescl*size(lyso_img,2)/size(lyso_img,1)/xscl];	yvals = ib(lysoinds)+sizescl*[0,1];
    for xc = 1:size(xvals,1)
      image(xvals(xc,:)-sizescl/2,yvals(xc,:)-sizescl/2,lyso_img,'AlphaData',lyso_alpha);
    end
  end

  if(~isempty(SUinds))
    xvals = allpos(SUinds)+[0,sizescl*size(SU_img,2)/size(SU_img,1)/xscl];	yvals = ib(SUinds)+sizescl*[0,1];
    for xc = 1:size(xvals,1)
      image(xvals(xc,:)-sizescl/2,yvals(xc,:)-sizescl/2,SU_img,'AlphaData',SU_alpha);
    end
  end

  if(~isempty(SFinds))
    xvals = allpos(SFinds)+[0,sizescl*size(SF_img,2)/size(SF_img,1)/xscl];	yvals = ib(SFinds)+sizescl*[0,1];
    for xc = 1:size(xvals,1)
      image(xvals(xc,:)-sizescl/2,yvals(xc,:)-sizescl/2,SF_img,'AlphaData',SF_alpha);
    end
  end

  if(~isempty(RUinds))
    xvals = allpos(RUinds)+[0,sizescl*size(RU_img,2)/size(RU_img,1)/xscl];	yvals = ib(RUinds)+sizescl*[0,1];
    for xc = 1:size(xvals,1)
      image(xvals(xc,:)-sizescl/2,yvals(xc,:)-sizescl/2,RU_img,'AlphaData',RU_alpha);
    end
  end

  if(~isempty(RFinds))
    xvals = allpos(RFinds)+[0,sizescl*size(RF_img,2)/size(RF_img,1)/xscl];	yvals = ib(RFinds)+sizescl*[0,1];
    for xc = 1:size(xvals,1)
      image(xvals(xc,:)-sizescl/2,yvals(xc,:)-sizescl/2,RF_img,'AlphaData',RF_alpha);
    end
  end
  plot([0,0]+sizescl/2,ylim,'k--')
  hold off

  plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28)
  ax = gca;
  ax.YTick = [];
  ax.TickLength = [0;0];
  pbaspect([pbscl,1,1])

  title(sprintf('time: %gs',round(tvals(snapind)*tau)))
  xlabel('distance from distal tip (\mum)')
  xlim([0,xmax])
  ylim([0,ymax])
end