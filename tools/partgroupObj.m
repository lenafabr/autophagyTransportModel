classdef partgroupObj < handle
	% object to store properties of groups of particles
	properties
		% number of particle types
		ntype = 2;
		
		% number of particles
		npart = 0;
		
		% list of particles
		particles = [];
		
		% snapshots
		snapshot = [];
		
	end
	methods
		function pgroup = partgroupObj
				pgroup.npart = [];
				pgroup.ntype = 2;
				pgroup.clearParticles;
		end
		
		function pgroup = clearParticles(pgroup,varargin)
			pgroup.particles = struct('time',[],'pos',[],'id',[],'prot',[],'acidity',[],'tdeg',[],...
											'xfuse',[],'tfuse',[],'tswitch',[],'nfuse',[],'type',[],'state',[],'intip',[]);
		end
	end
end