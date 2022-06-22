## Functions to run energy evaluation in LAMMPS.jl
using LAMMPS

function startN_LAMMPS_instances(N::Int)
	LMPvect=[ LMP(["-screen","none","-sf", "omp", "-pk", "omp 6"]) for i=1:N]
	#LMPvect=[ LMP(["-sf", "omp", "-pk", "omp 4"]) for i=1:N]
	return LMPvect
end



function MinimizeCoords(Full_LMPvect,OLatticeSites,HFLatticeSites,Temp,MD_timestep,SimDim,MinSteps,MDsteps)
	command.(LMPvect,"clear")
	OxygenAtoms=OLatticeSites[OLatticeSites[:,7].==1,:];
	

	LMP_v=Full_LMPvect[1];

	FinBaseCoords=[[Inf Inf Inf Inf] for ind=1:size(HFLatticeSites,1)];
	FinOxyCoords=[[Inf Inf Inf Inf] for ind=1:size(OxygenAtoms,1)];
	newOLatticeSites=[];
	newHFLatticeSites=[];
	#Basic LMP setup
	try
		command(LMP_v,"units  metal");
		command(LMP_v,"atom_style  charge");
		command(LMP_v,"dimension  3");
		command(LMP_v,"boundary  p p m");
		command(LMP_v,"processors * * *");
		#Setup Box
		command(LMP_v,"region SimulationDomain block 0.0 " * string(SimDim[1]) * " 0.0 " * string(SimDim[2]) * " 0.0 " * string(SimDim[3]));
		command(LMP_v,"create_box 4 SimulationDomain");
		
		for jnd=1:size(HFLatticeSites,1)
			command(LMP_v,"create_atoms " * string(floor(Int,HFLatticeSites[jnd,2])) * " single " * string(HFLatticeSites[jnd,4]) * " " * string(HFLatticeSites[jnd,5]) * " " * string(HFLatticeSites[jnd,6]) )
		end
		for knd=1:size(OxygenAtoms,1)
			command(LMP_v,"create_atoms 4 single " * string(OxygenAtoms[knd,4]) * " " * string(OxygenAtoms[knd,5]) * " " * string(OxygenAtoms[knd,6]));
		end
		#Continue Setup
		command(LMP_v,"group Hf-fixed type 1");
		command(LMP_v,"group Hf-temp type 2");
		command(LMP_v,"group Hf-nve type 3");
		command(LMP_v,"group Oxygen 	type 4");

		command(LMP_v,"group Hf-group 	type 1 2 3");
		
		
		command(LMP_v,"mass 1 178.4900");
		command(LMP_v,"mass 2 178.4900");
		command(LMP_v,"mass 3 178.4900");
		command(LMP_v,"mass 4 15.9990");

		#Setup Potential
		command(LMP_v,"pair_style  comb");
		command(LMP_v,"pair_coeff  * * ffield.comb Hf Hf Hf O ");
		command(LMP_v,"neighbor  2.0 bin");
		command(LMP_v,"neigh_modify  every 1 delay 0 check no");

		command(LMP_v,"fix  CombQEQ    all  qeq/comb  1  0.001");
		#println("Set QEQcomb: minimize");

		
		command(LMP_v,"timestep  " *  string(MD_timestep));
		
		#Define Dynamics
		command(LMP_v,"velocity  all create " * string(Temp) * " " * string(abs.(rand(Int16,1))[1]) * " dist gaussian");
		
		command(LMP_v,"fix  101  Hf-fixed  move linear  0.0  0.0  0.0  units box ");

		command(LMP_v,"fix  102 Hf-temp   nvt temp " * string(Temp) * " " * string(Temp) * " "  * string(50.0*MD_timestep));
		command(LMP_v,"fix  103 Hf-nve   nve");
		command(LMP_v,"fix  104 Oxygen   nve");
		
		#command(LMP_v,"delete_atoms overlap 0.5 all all")
		
		command(LMP_v,"compute  AtomPE all pe/atom ");
		command(LMP_v,"compute  AtomXYZ all property/atom x y z ");
		command(LMP_v,"compute  MyTestPE all reduce sum c_AtomPE");
		command(LMP_v,"compute  MyTestLoc all reduce sum c_AtomXYZ[1] c_AtomXYZ[2] c_AtomXYZ[3]");
		command(LMP_v,"thermo 	20");
		command(LMP_v,"thermo_style 	custom step temp press etotal atoms c_MyTestPE spcpu");
		command(LMP_v,"min_style sd");
		command(LMP_v,"min_modify dmax 0.01");
		command(LMP_v,"minimize	1.0e-10  1.0e-10  " * string(MinSteps) * " "  * string(MinSteps));
		command(LMP_v,"run " * string(MDsteps));
		
		FinCoords=hcat(extract_atom(LMP_v,"type"),convert(Matrix,extract_atom(LMP_v,"x")'))

		OxyKey=FinCoords[:,1].==4;
		BaseKey=.!OxyKey;
		newOLatticeSites=convert(Matrix{Any},hcat(1:sum(OxyKey),FinCoords[OxyKey,1],zeros(sum(OxyKey),1),FinCoords[OxyKey,2:4],ones(sum(OxyKey),1)));
		newHFLatticeSites=convert(Matrix{Any},hcat(1:sum(BaseKey),FinCoords[BaseKey,1],zeros(sum(BaseKey),1),FinCoords[BaseKey,2:4],ones(sum(BaseKey),1)));
	catch
		newOLatticeSites=OLatticeSites;
		newHFLatticeSites=HFLatticeSites;
		println("Minimize Sim Failed");
	end
	
	return newOLatticeSites,newHFLatticeSites
end

function MinimizeListedAtoms(Full_LMPvect,OLatticeSites,HFLatticeSites,MinOxyVec,MinHfVec,Temp,MD_timestep,SimDim,MinSteps,MDsteps)
	command.(LMPvect,"clear")
	OxygenAtoms=OLatticeSites[OLatticeSites[:,7].==1,:];
	
	if !isempty(MinOxyVec)
		O_key = in(MinOxyVec).(OLatticeSites[:,1])
	else
		O_key = zeros(Bool,size(OLatticeSites,1));
	end
	if !isempty(MinHfVec)	
		Hf_key= in(MinHfVec).(HFLatticeSites[:,1])
	else
		Hf_key = zeros(Bool,size(HFLatticeSites,1))
	end
	
	OLat_min=OLatticeSites[O_key,:];
	HfLat_min=HFLatticeSites[Hf_key,:];
	OLat_freeze=OLatticeSites[.!O_key,:];
	HfLat_freeze=HFLatticeSites[.!Hf_key,:];
	
	LMP_v=Full_LMPvect[1];

	FinBaseCoords=[[Inf Inf Inf Inf] for ind=1:size(HFLatticeSites,1)];
	FinOxyCoords=[[Inf Inf Inf Inf] for ind=1:size(OxygenAtoms,1)];
	newOLatticeSites=[];
	newHFLatticeSites=[];
	#Basic LMP setup
	try
		command(LMP_v,"units  metal");
		command(LMP_v,"atom_style  charge");
		command(LMP_v,"dimension  3");
		command(LMP_v,"boundary  p p m");
		command(LMP_v,"processors * * *");
		#Setup Box
		command(LMP_v,"region SimulationDomain block 0.0 " * string(SimDim[1]) * " 0.0 " * string(SimDim[2]) * " 0.0 " * string(SimDim[3]));
		command(LMP_v,"create_box 8 SimulationDomain");
		
		for jnd=1:size(HfLat_freeze,1)
			command(LMP_v,"create_atoms " * string(floor(Int,HfLat_freeze[jnd,2])) * " single " * string(HfLat_freeze[jnd,4]) * " " * string(HfLat_freeze[jnd,5]) * " " * string(HfLat_freeze[jnd,6]) )
		end
		for knd=1:size(OLat_freeze,1)
			command(LMP_v,"create_atoms 4 single " * string(OLat_freeze[knd,4]) * " " * string(OLat_freeze[knd,5]) * " " * string(OLat_freeze[knd,6]));
		end
		for jnd=1:size(HfLat_min,1)
			command(LMP_v,"create_atoms " * string(floor(Int,HfLat_min[jnd,2].+4)) * " single " * string(HfLat_min[jnd,4]) * " " * string(HfLat_min[jnd,5]) * " " * string(HfLat_min[jnd,6]) )
		end
		for knd=1:size(OLat_min,1)
			command(LMP_v,"create_atoms 8 single " * string(OLat_min[knd,4]) * " " * string(OLat_min[knd,5]) * " " * string(OLat_min[knd,6]));
		end		
		
		#Continue Setup
		command(LMP_v,"group Hf-fixed type 1 2 3 5");
		command(LMP_v,"group Hf-temp type 6");
		command(LMP_v,"group Hf-nve type 7");
		command(LMP_v,"group Oxygen 	type 8");
		command(LMP_v,"group O-fixed type 4");

		command(LMP_v,"group Hf-group 	type 1 2 3 5 6 7");
		command(LMP_v,"group all-fixed 	type 1 2 3 4 5");
		
		command(LMP_v,"mass 1 178.4900");
		command(LMP_v,"mass 2 178.4900");
		command(LMP_v,"mass 3 178.4900");
		command(LMP_v,"mass 4 15.9990");
		command(LMP_v,"mass 5 178.4900");
		command(LMP_v,"mass 6 178.4900");
		command(LMP_v,"mass 7 178.4900");
		command(LMP_v,"mass 8 15.9990");
		
		#Setup Potential
		command(LMP_v,"pair_style  comb");
		command(LMP_v,"pair_coeff  * * ffield.comb Hf Hf Hf O Hf Hf Hf O");
		command(LMP_v,"neighbor  2.0 bin");
		command(LMP_v,"neigh_modify  every 1 delay 0 check no");

		command(LMP_v,"fix  CombQEQ    all  qeq/comb  1  0.001");
		#println("Set QEQcomb: minimize");

		
		command(LMP_v,"timestep  " *  string(MD_timestep));
		
		#Define Dynamics
		command(LMP_v,"velocity  all create " * string(Temp) * " " * string(abs.(rand(Int16,1))[1]) * " dist gaussian");
		
		command(LMP_v,"fix  101  all-fixed  move linear  0.0  0.0  0.0  units box ");

		command(LMP_v,"fix  102 Hf-temp   nvt temp " * string(Temp) * " " * string(Temp) * " "  * string(50.0*MD_timestep));
		command(LMP_v,"fix  103 Hf-nve   nve");
		command(LMP_v,"fix  104 Oxygen   nve");
		
		#command(LMP_v,"delete_atoms overlap 0.5 all all")
		
		command(LMP_v,"compute  AtomPE all pe/atom ");
		command(LMP_v,"compute  AtomXYZ all property/atom x y z ");
		command(LMP_v,"compute  MyTestPE all reduce sum c_AtomPE");
		command(LMP_v,"compute  MyTestLoc all reduce sum c_AtomXYZ[1] c_AtomXYZ[2] c_AtomXYZ[3]");
		command(LMP_v,"thermo 	20");
		command(LMP_v,"thermo_style 	custom step temp press etotal atoms c_MyTestPE spcpu");
		command(LMP_v,"min_style sd");
		command(LMP_v,"min_modify dmax 0.01");
		command(LMP_v,"minimize	1.0e-10  1.0e-10  " * string(MinSteps) * " "  * string(MinSteps));
		command(LMP_v,"run " * string(MDsteps));
		
		FinCoords=hcat(extract_atom(LMP_v,"type"),convert(Matrix,extract_atom(LMP_v,"x")'))

		OxyKey=(FinCoords[:,1].==4) .|| (FinCoords[:,1].==8);
		BaseKey=.!OxyKey;
		newOLatticeSites=convert(Matrix{Any},hcat(1:sum(OxyKey),FinCoords[OxyKey,1],zeros(sum(OxyKey),1),FinCoords[OxyKey,2:4],ones(sum(OxyKey),1)));
		newOLatticeSites[:,2].=4;
		newHFLatticeSites=convert(Matrix{Any},hcat(1:sum(BaseKey),FinCoords[BaseKey,1],zeros(sum(BaseKey),1),FinCoords[BaseKey,2:4],ones(sum(BaseKey),1)));
		KeyFrozen=newHFLatticeSites[:,2].>4;
		newHFLatticeSites[KeyFrozen,2]=newHFLatticeSites[KeyFrozen,2].-4;
		
	catch
		newOLatticeSites=OLatticeSites;
		newHFLatticeSites=HFLatticeSites;
		println("Minimize Sim Failed");
	end
	
	return newOLatticeSites,newHFLatticeSites
end



function EnergyEvalSim(Full_LMPvect,OLatticeSites,HFLatticeSites,PossibleNeighbors,Temp,MD_timestep,SimDim)
	command.(LMPvect,"clear")
	OxygenAtoms=OLatticeSites[OLatticeSites[:,7].==1,:];
	NumSites=size(PossibleNeighbors,1);
	LMP_v=Full_LMPvect[1:NumSites];
	EnergyVector=[Inf for ind=1:NumSites]

	#Basic LMP setup
	#Threads.@threads 
	for ind=1:NumSites
		
		try
			command(LMP_v[ind],"units  metal");
			command(LMP_v[ind],"atom_style  charge");
			command(LMP_v[ind],"dimension  3");
			command(LMP_v[ind],"boundary  p p m");
			command(LMP_v[ind],"processors * * *");
			#Setup Box
			command(LMP_v[ind],"region SimulationDomain block 0.0 " * string(SimDim[1]) * " 0.0 " * string(SimDim[2]) * " 0.0 " * string(SimDim[3]));
			command(LMP_v[ind],"create_box 5 SimulationDomain");
			

			command(LMP_v[ind],"create_atoms 5 single " * string(PossibleNeighbors[ind,1]) * " " * string(PossibleNeighbors[ind,2]) * " " * string(PossibleNeighbors[ind,3]));

			for jnd=1:size(HFLatticeSites,1)
				command(LMP_v[ind],"create_atoms " * string(floor(Int,HFLatticeSites[jnd,2])) * " single " * string(HFLatticeSites[jnd,4]) * " " * string(HFLatticeSites[jnd,5]) * " " * string(HFLatticeSites[jnd,6]) )
			end
			for knd=1:size(OxygenAtoms,1)
				command(LMP_v[ind],"create_atoms 4 single " * string(OxygenAtoms[knd,4]) * " " * string(OxygenAtoms[knd,5]) * " " * string(OxygenAtoms[knd,6]));
			end
			#Continue Setup
			command(LMP_v[ind],"group Hf-fixed type 1");
			command(LMP_v[ind],"group Hf-temp type 2");
			command(LMP_v[ind],"group Hf-nve type 3");
			command(LMP_v[ind],"group Oxygen 	type 4");
			command(LMP_v[ind],"group OxygenTest 	type 5");
			command(LMP_v[ind],"group Freeze 	type 1 2 3 4");
			command(LMP_v[ind],"mass 1 178.4900");
			command(LMP_v[ind],"mass 2 178.4900");
			command(LMP_v[ind],"mass 3 178.4900");
			command(LMP_v[ind],"mass 4 15.9990");
			command(LMP_v[ind],"mass 5 15.9990");
			#Setup Potential
			command(LMP_v[ind],"pair_style  comb");
			command(LMP_v[ind],"pair_coeff  * * ffield.comb Hf Hf Hf O O ");
			command(LMP_v[ind],"neighbor  2.0 bin");
			command(LMP_v[ind],"neigh_modify  every 1 delay 0 check no");
			
			command(LMP_v[ind],"fix  CombQEQ    all  qeq/comb  1  0.001");
			
			println("Set QEQcomb: sim" * string(ind));
			
			command(LMP_v[ind],"timestep  " *  string(MD_timestep));
			#println("Set timestep");
			#Define Dynamics
			#command(LMP_v[ind],"delete_atoms overlap 0.5 all all")
			command(LMP_v[ind],"velocity  all create " * string(Temp) * " " * string(abs.(rand(Int16,1))[1]) * " dist gaussian");
			#println("Set velocity");
			command(LMP_v[ind],"fix  FreezeAtoms  Freeze  move linear  0.0  0.0  0.0  units box ");
			command(LMP_v[ind],"fix  102 OxygenTest   nve/limit 0.05");
			#println("Set fix");
			command(LMP_v[ind],"compute  AtomPE OxygenTest pe/atom ");
			command(LMP_v[ind],"compute  AtomXYZ OxygenTest property/atom x y z ");
			command(LMP_v[ind],"compute  MyTestPE OxygenTest reduce sum c_AtomPE");
			command(LMP_v[ind],"compute  MyTestLoc OxygenTest reduce sum c_AtomXYZ[1] c_AtomXYZ[2] c_AtomXYZ[3]");
			#println("Set compute");
			command(LMP_v[ind],"thermo 	1");
			command(LMP_v[ind],"thermo_style 	custom step temp press etotal atoms c_MyTestPE spcpu");
			command(LMP_v[ind],"min_style sd");
			command(LMP_v[ind],"min_modify dmax 0.1");
			#println("MINIMIZE Energy Eval")
			command(LMP_v[ind],"minimize	1.0e-10  1.0e-10  15  15");
			#println("MD Energy Eval")
			command(LMP_v[ind],"run 10");
			
			atomIDkey=extract_atom(LMP_v[ind],"type").==5
			EnergyVector[ind]=extract_compute(LMP_v[ind],"AtomPE",LAMMPS.API.LMP_STYLE_ATOM,LAMMPS.API.LMP_TYPE_VECTOR)[atomIDkey][1]
		catch
			#println("Bad LAMMPS Sim Caught")
			EnergyVector[ind]=Inf;
		
			
		end
	end
	#display(FinCoords)
	#display(EnergyVector)
	command.(LMPvect,"clear")
	
	#command.(LMP_v,"");
	#command.(LMPvect,"clear")
	
	return EnergyVector
end


