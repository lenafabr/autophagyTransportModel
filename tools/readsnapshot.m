function [grouplist,tvals,domlen,ntrials] = readsnapshot(filename,options)
% READSNAPSHOT reads out particle states and positions from filename
% *.snap.out file must be created using partdynamics1D.exe
% output format for each timestep:
% Info line: TIME NTRIALS DOMLEN
%						|Ntrial blocks: TRIALNUMBER NPART
%						|							|Npart blocks: pc, ID, type, nprot, nfuse, state, intip
%						|							|							pos, prot, acidity, degtime, xfuse,	tfuse, tswitch

opt = struct();
% output progress after how many steps
opt.progstep = 0;
% only read a specific particle type
opt.onlyread = 0;
% maximum number of particles to read
opt.maxpart = inf;
% only read last snapshot?
opt.readlastsnap = 0;
% get snapshots for movies?
opt.getmoviesnaps = 0;

if(exist('options','var')); opt = copyStruct(options,opt); end

fid = fopen(filename,'r');
tvals = [];
snapct = 0;

while(~feof(fid))
	initdata = textscan(fid,'%f %f %f',1); % time ntrials domlen
	if(isempty(initdata{1})); return; end
	snapct = snapct+1;
	if(~mod(snapct,opt.progstep));  disp(snapct); end
	curtime = initdata{1};
	tvals = cat(1,tvals,curtime);
	if(~exist('ntrials','var')); ntrials = initdata{2}; end
	if(~exist('domlen','var')); domlen = initdata{3}; end
	if(~exist('grouplist','var')); grouplist(ntrials,1) = partgroupObj; end

	% start reading all trials for a single time step
	trc = 1;
	while(trc<=ntrials)
		trinfo = textscan(fid,'%f %f',1); % trialnumber npart
		tc = trinfo{1};
		npart = trinfo{2};
		
		grouplist(tc).npart = cat(1,grouplist(tc).npart,npart);
		
		% read all particles for a given trial
		pdata = textscan(fid,'%f %f %f %f %f %f %f %f',2*npart); %pc ID type nprot nfuse state intip
																									  %pos prot acidity degtime xfuse tfuse tswitch
		if(npart<=0)
			trc = trc+1;
			continue;
		end
		pidvals = pdata{2}(1:2:end);
		typevals = pdata{3}(1:2:end);
		nfusevals = pdata{5}(1:2:end);
		statevals = pdata{6}(1:2:end);
    intipvals = pdata{7}(1:2:end);
		
		posvals = pdata{1}(2:2:end);
		protvals = [pdata{2}(2:2:end),pdata{3}(2:2:end)];
		acidityvals = pdata{4}(2:2:end);
		tdegvals = pdata{5}(2:2:end);
		xfusevals = pdata{6}(2:2:end);
		tfusevals = pdata{7}(2:2:end);
		tswitchvals = pdata{8}(2:2:end);
		
		if(opt.getmoviesnaps)
			grouplist(tc).snapshot(snapct).time = curtime;
			grouplist(tc).snapshot(snapct).id = pidvals;
			grouplist(tc).snapshot(snapct).type = typevals;
			grouplist(tc).snapshot(snapct).nfuse = nfusevals;
			grouplist(tc).snapshot(snapct).state = statevals;
      grouplist(tc).snapshot(snapct).intip = intipvals;
			grouplist(tc).snapshot(snapct).pos = posvals;
			grouplist(tc).snapshot(snapct).prot = protvals;
			grouplist(tc).snapshot(snapct).acidity = acidityvals;
			grouplist(tc).snapshot(snapct).tdeg = tdegvals;
			grouplist(tc).snapshot(snapct).xfuse = xfusevals;
			grouplist(tc).snapshot(snapct).tfuse = tfusevals;
			grouplist(tc).snapshot(snapct).tswitch = tswitchvals;
		end
		
		if(opt.readlastsnap)
			grouplist(tc).npart = npart;
			grouplist(tc).particles = [];
			for pc = 1:length(pidvals)
				grouplist(tc).particles(pc).id = pidvals(pc);
				grouplist(tc).particles(pc).type = typevals(pc);
				grouplist(tc).particles(pc).time = initdata{1};
				grouplist(tc).particles(pc).pos = posvals(pc);
				grouplist(tc).particles(pc).prot = protvals(pc,:);
				grouplist(tc).particles(pc).nfuse = nfusevals(pc);
				grouplist(tc).particles(pc).acidity = acidityvals(pc);
				grouplist(tc).particles(pc).tdeg = tdegvals(pc);
				grouplist(tc).particles(pc).xfuse = xfusevals(pc);
				grouplist(tc).particles(pc).tfuse = tfusevals(pc);
				grouplist(tc).particles(pc).tswitch = tswitchvals(pc);				
				grouplist(tc).particles(pc).state = statevals(pc);
        grouplist(tc).particles(pc).intip = intipvals(pc);
			end
		else
			for pc = 1:length(pidvals)
				pid = pidvals(pc);
				grouplist(tc).particles(pid).id = pid;
				grouplist(tc).particles(pid).type = typevals(pc);
				grouplist(tc).particles(pid).time = [grouplist(tc).particles(pid).time;initdata{1}];
				grouplist(tc).particles(pid).pos = [grouplist(tc).particles(pid).pos;posvals(pc)];
				grouplist(tc).particles(pid).prot = [grouplist(tc).particles(pid).prot;protvals(pc,:)];
				grouplist(tc).particles(pid).nfuse = nfusevals(pc);
				grouplist(tc).particles(pid).acidity = [grouplist(tc).particles(pid).acidity;acidityvals(pc)];
				grouplist(tc).particles(pid).tdeg = [grouplist(tc).particles(pid).tdeg;tdegvals(pc)];
				grouplist(tc).particles(pid).xfuse = [grouplist(tc).particles(pid).xfuse;xfusevals(pc)];
				grouplist(tc).particles(pid).tfuse = [grouplist(tc).particles(pid).tfuse;tfusevals(pc)];
				grouplist(tc).particles(pid).tswitch = [grouplist(tc).particles(pid).tswitch;tswitchvals(pc)];
				grouplist(tc).particles(pid).state = [grouplist(tc).particles(pid).state;statevals(pc)];
        grouplist(tc).particles(pid).intip = [grouplist(tc).particles(pid).intip;intipvals(pc)];
			end
		end
		
		trc = trc+1;
	end
end

fclose(fid);
end

