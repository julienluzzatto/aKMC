#!/usr/bin/env julia

function writeMinimization(HfAtoms,OxAtoms,BoxSizes,file="Atoms.dat")

	open("in.min", "w") do io
		println(io, "clear");

		# Basic setup
		println(io, "units metal");
		#println(io, "atom_style atomic");
		println(io, "atom_style charge");
		println(io, "atom_modify map array sort 0 0.0");
		println(io, "dimension 3");
		println(io, "boundary p p m");
		println(io, "processors * * *");

		# Create box
		println(io, "region simdim block 0.0 " * string(BoxSizes[1]) * " 0.0 " * string(BoxSizes[2]) * " 0.0 " * string(BoxSizes[3]));
		OxPositions = pyconvert(Matrix{Float64},OxAtoms.positions);
		if size(OxPositions,1) > 0
			println(io, "create_box 2 simdim");
		else
			println(io, "create_box 1 simdim");
		end

		# Import configuration
		println(io, "read_data " * file * " add append");

		# Define regular groups
		println(io, "group Hf type 1");
		if size(OxPositions,1) > 0
			println(io, "group Ox type 2");
		end

		# Define properties
		println(io, "mass 1 178.4900");
		if size(OxPositions,1) > 0
			println(io, "mass 2 15.9990");
		end

		# Define potential
		#println(io, "pair_style pod");
		println(io, "pair_style comb");
		if size(OxPositions,1) > 0
			#println(io, "pair_coeff * * pod.txt coefficients.txt Hf O");
			println(io, "pair_coeff * * ffield.comb Hf O");
		else
			#println(io, "pair_coeff * * pod.txt coefficients.txt Hf");
			println(io, "pair_coeff * * ffield.comb Hf");
		end
		println(io, "fix CombQEQ all qeq/comb 1 0.001");
		println(io, "neighbor 1.0 bin");
		println(io, "neigh_modify every 1 delay 0 check yes");

		# Dump new coordinates
		println(io, "dump 1 Hf custom 100000 dump.Hf id type x y z");
		if size(OxPositions,1) > 0
			println(io, "dump 2 Ox custom 100000 dump.Ox id type x y z");
		end

		# Minimize coordinates
		println(io, "fix 1 all box/relax x 0.0 y 0.0");
		println(io, "min_style cg");
		println(io, "min_modify dmax 1.0e-2 line quadratic");
		println(io, "minimize 1e-10 1e-10 10000 10000");
	end
end

function writeNEB(HfAtoms,OxAtoms,BoxSizes,Target,file="Atoms.dat")
#function NEB_write_input(OLatticeSites,HFLatticeSites,Temp,MD_timestep,SimDim,MinSteps,MDsteps,target,coords)

	open("in.neb", "w") do io
		println(io, "clear");

		# Basic setup
		println(io, "units metal");
		#println(io, "atom_style atomic");
		println(io, "atom_style charge");
		println(io, "atom_modify map array sort 0 0.0");
		println(io, "dimension 3");
		println(io, "boundary p p m");
		println(io, "processors * * *");

		# Create box
		println(io, "region simdim block 0.0 " * string(BoxSizes[1]) * " 0.0 " * string(BoxSizes[2]) * " 0.0 " * string(BoxSizes[3]));
		OxPositions = pyconvert(Matrix{Float64},OxAtoms.positions);
		println(io, "create_box 2 simdim");

		# Import configuration
		println(io, "read_data " * file * " add append");

		# Define regular groups
		println(io, "group Hf type 1");
		println(io, "group Ox type 2");

		# Define properties
		println(io, "mass 1 178.4900");
		println(io, "mass 2 15.9990");

		# Define potential
		#println(io, "pair_style pod");
		println(io, "pair_style comb");
		#println(io, "pair_coeff * * pod.txt coefficients.txt Hf O");
		println(io, "pair_coeff * * ffield.comb Hf O");
		println(io, "fix CombQEQ all qeq/comb 1 0.001");
		println(io, "neighbor 0.3 bin");
		println(io, "neigh_modify every 1 delay 5");

		# Define nudged elastic band groups
		println(io, "group atomNEB id " * string(Target));
		println(io, "group noneNEB subtract all atomNEB");

		# Minimize before NEB
		println(io, "min_style cg");
		println(io, "minimize 1e-5 1e-5 100 100");

		# Characterize nudged elastic band
		println(io, "fix 1 atomNEB neb 1.0");
		println(io, "fix 2 noneNEB setforce 0.0 0.0 0.0");
		println(io, "thermo 1");

		# Dump new coordinates
		println(io, "dump 1 atomNEB custom 10 dump.neb id type x y z")

		# Run nudged elastic band calculation
		println(io, "min_style quickmin");
		println(io, "neb 1e-3 1e-4 100 100 10 final coords.dat");
	end
end
