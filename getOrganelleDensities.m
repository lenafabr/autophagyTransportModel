% GETORGANELLEDENSITIES Calculates the mean-field densities of autophagosomes (AVs) and
% lysosomes involved in neuronal autophagy.
%% set up options for calculation
options = struct();

options.nx = 1000; % number of points on the calculation grid

options.solvePfuse = true;
options.dotipfusion = true; % fuse particles in the distal tip
options.onedone = false; % allow AVs to only fuse once

numBranch = 5; % number of branches

xvals = linspace(0,1,options.nx);
dx = xvals(2)-xvals(1);
%% set up parameter structure
params = struct();

% domain length and time scales
params.Lreal = 1055; % length of domain in real units (um). Sets length scale
params.vpreal = 0.75; % AV velocity in real units (um/s). Sets time scale
params.L = 1; % scaled domain length
params.vp = 1; % scaled AV velocity

% organelle sizes
params.rp = 0.4/params.Lreal; % AV radius
params.ry = 0.1/params.Lreal; % lysosome radius
params.ell = 2*(params.rp+params.ry); % minimum contact distance

% AV production params
params.kpp0 = 1; % AV production rate in main axon tip
params.kpp = params.kpp0*ones(1,numBranch); % AV production rate in branch tips

% params for bidirectional AV motion
params.vpb = 1; % AV velocity in bidirectional state
lam = 2; % AV run length in bidirectional state (um)
run_length = lam/params.Lreal; % scaled run length
fstat = 0.4603; % fraction of stationary AVs (obtained from exp data)
params.khb = params.vpb/run_length; % rate of halting in the bidirectional state
params.kwb = params.khb*(1-fstat)/fstat; % rate of walking in the bidirectional state
params.pretro = 0.5; % probability of walking retrograde upon restarting

% params for persistently retrograde AV motion
params.ks = 1; % rate of AVs switching to persistent retrograde motion
params.kwr = 1; % rate of walking in the retrograde state
haltratio = 0.6578; % ratio of halting rate to walking rate for retrograde AVs
params.khr = haltratio*params.kwr; % rate of halting in the retrograde state

% parameters for lysosome production and motion
params.kpy = 120; % lysosome production rate
% 40% retro @ 1.44um/s, 30% stationary, 30% antero @ 2um/s
% source: Boecker et al. 2020
params.vya = 0.3*2/(0.3+0.15)/params.vpreal; % anterograde lysosome velocity
params.vyr = 0.4*1.5/(0.4+0.15)/params.vpreal; % retrograde lysosome velocity

params.khy = 0; %lysosome halting rate
params.kwy = 0; %lysosome walking rate

% parameters for distal tip
params.kye = 1; % rate of lysosome return from distal tip
params.kpe = inf; % rate of AV return from distal tip

% parameters for lysosome decay
params.kda = 0; params.kdr = 0; params.kds = 0;

% parameters for fusion probability
params.pf = 0.005;
% first letter denotes AV direction, second letter denotes lysosome direction
% a: anterograde, r: retrograde, s: stationary
params.pfra = params.pf; params.pfrr = params.pf; params.pfrs = params.pf;
params.pfaa = params.pf; params.pfar = params.pf; params.pfas = params.pf;
params.pfsa = params.pf; params.pfsr = params.pf;
params.pft = params.pf; % fusion probability within tip

% options and parameters for axonal branching
options.pfsplit = false; %same fusion probability in all branches
params.mainrad = 1*ones(1,numBranch+1); % radii of main axon segments
params.collrad = 1*ones(1,numBranch+1); % radii of axon collaterals/branches
params.branchpoints = linspace(250,params.Lreal-250,numBranch)/params.Lreal; % location of branch points along the main axon
params.branchlen = 164/params.Lreal*ones(size(params.branchpoints)); % branch lengths

% split organelles according to number of branch tips served
params.psplit = [];
for c = 1:numBranch
  params.psplit(c) = c/(c+1);
end
%% obtain density structure
if(options.onedone)
  [dens,success] = onedone_branches_solver(params,options);
else
  [dens,success] = unlim_branches_solver(params,options);
end
if(success); fprintf('Densities calculated successfully.\n'); end
