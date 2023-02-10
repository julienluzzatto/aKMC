#!/usr/bin/env julia

#global
using ASEconvert
using LinearAlgebra
using Interpolations
using UnicodePlots
using DelimitedFiles
using CMAEvolutionStrategy
using Glob
using StatsBase
using QHull
using MiniQhull
using GeometricalPredicates
using LoopVectorization
using AtomsBase

command1 = `mpirun -np 32 /home/gridsan/jluzzatto/lammps/build/lmp -screen none -in in.min`
command2 = `mpirun -np 32 /home/gridsan/jluzzatto/lammps/build/lmp -partition 8x4 -screen none -in in.neb`;

function makeExistingConfiguration(Dimensions,Structure=0,LatCst=3.5416)
	if Structure == 0 # BCC crystal structure
		HfAtoms = ase.build.bcc100("Hf",Dimensions,a=LatCst);
        OxAtoms = ase.Atoms();
	else # HCP crystal Structure
		HfAtoms = ase.build.hcp0001("Hf",Dimensions,orthogonal=True)
        OxAtoms = ase.Atoms();
	end
	return HfAtoms,OxAtoms
end

function loadExistingConfiguration()
	HfAtoms = load_system("HfAtoms.xyz");
	OxAtoms = load_system("OxAtoms.xyz");
	return HfAtoms,OxAtoms
end

function getDimensions(Dimensions,Structure=0,LatCst=3.5416)
	if Structure == 0 # BCC crystal structure
		UnitCell = [LatCst,LatCst,LatCst];
		BoxSizes = Dimensions .* UnitCell;
	else # HCP crystal Structure
		UnitCell = [LatCst,LatCst,LatCst];
		BoxSizes = Dimensions .* UnitCell; # inaccurate for HCP
	end
	return UnitCell,BoxSizes
end

function runMinimization(HfAtoms,OxAtoms,BoxSizes)
	ase.io.write("Atoms.dat",HfAtoms+OxAtoms,format="lammps-data",atom_style="charge");
	writeMinimization(HfAtoms,OxAtoms,BoxSizes);
	command = command1;
	read(command)
	HfAtoms,OxAtoms,BoxSizes = MinimizeCoords(HfAtoms,OxAtoms,BoxSizes);
end

function MinimizeCoords(HfAtoms,OxAtoms,BoxSizes)
	Haf = ase.io.read("dump.Hf",format="lammps-dump-text");
	Haf.set_chemical_symbols(["Hf" for at in Haf]);
	Haf.wrap(pbc=[1,1,0]);

	OxPositions = pyconvert(Matrix{Float64},OxAtoms.positions);
	if size(OxPositions,1) > 0
		Oxy = ase.io.read("dump.Ox",format="lammps-dump-text");
		Oxy.set_chemical_symbols(["O" for at in Oxy]);
		Oxy.wrap(pbc=[1,1,0]);
	else
		Oxy = ase.Atoms();
	end

	Box = pyconvert(Array{Float64},[Haf.cell[0,0],Haf.cell[1,1],Haf.cell[2,2]]);
	println("Minimization Achieved");
	return Haf,Oxy,Box
end

function GenerateOxygenTrialPoints(Box,Unit,HfLattice,MinSpacing)
    xSize,ySize,zSize = Box;
	BaseHaf = pyconvert(Matrix{Any},HfLattice.positions)

    # Create ghost atoms to account for periodic boundary conditions
    extBaseLattice = BaseHaf[:,1:3];
    BLxlo = extBaseLattice[extBaseLattice[:,1].<4*Unit[1],:].+[xSize,0,0]';
    BLxhi = extBaseLattice[extBaseLattice[:,1].>xSize-(4*Unit[1]),:].-[xSize,0,0]';
    extBaseLattice = [extBaseLattice;BLxlo;BLxhi];
    BLylo = extBaseLattice[extBaseLattice[:,2].<4*Unit[2],:].+[0,ySize,0]';
    BLyhi = extBaseLattice[extBaseLattice[:,2].>ySize-(4*Unit[2]),:].-[0,ySize,0]';
    extBaseLattice = [extBaseLattice;BLylo;BLyhi];
    BLzlo = extBaseLattice[(extBaseLattice[:,3].<3*Unit[3]) .& (extBaseLattice[:,3].>0),:].*[1,1,-1]';
    extBaseLattice = [extBaseLattice;BLzlo];
    numBasePoints = size(extBaseLattice,1);

    #Account for surface roughness
    #ID surface atoms
    function nlargest(v, n; rev=true)
        result = falses(size(v))
        result[partialsortperm(v, 1:n; rev=rev)] .= true
        return result
    end
	RepZ = round(zSize/Unit[3]);
    KEY_upperExtBL = nlargest(extBaseLattice[:,3],minimum([round(Int,numBasePoints*2/RepZ),numBasePoints]));

    Index_upperExtBL = collect(1:1:numBasePoints)
    Ind_surfaceAtoms = Index_upperExtBL[KEY_upperExtBL];
    surfaceAtoms = extBaseLattice[KEY_upperExtBL,:];
    keepAtom = zeros(Int,size(Ind_surfaceAtoms,1));

    for ind = 1:size(Ind_surfaceAtoms,1)
        locPos = surfaceAtoms[ind,:];
        searchSpacing = maximum(Unit[1:2])*3/8;
        xdist = surfaceAtoms[:,1].-locPos[1];
        ydist = surfaceAtoms[:,2].-locPos[2];
        neighKEY = (abs.(xdist).<searchSpacing).&(abs.(ydist).<searchSpacing);
        neighZheights = surfaceAtoms[neighKEY,3];

        isheighest = locPos[3]>=maximum(neighZheights); #loc is highest? true false
        keepAtom[ind] = isheighest;
    end
    surfaceAtoms = surfaceAtoms[convert.(Bool,keepAtom),:];
    Ind_surfaceAtoms = Ind_surfaceAtoms[convert.(Bool,keepAtom)];
    surfDelaunayPoly = delaunay(2,size(surfaceAtoms,1),vec(surfaceAtoms[:,1:2]'))

    #Calcluate Regular Delaunay Tess
    delaunayPoly = delaunay(3,numBasePoints,vec(extBaseLattice[:,1:3]'));
    midMatrix = zeros(size(delaunayPoly,2)*6,3);
    thirdMatrix = zeros(size(delaunayPoly,2)*4,3);
    fourthMatrix = zeros(size(delaunayPoly,2)*1,3);

    #Calculate 2,3, and 4 fold midpoints
    for ind = 1:size(delaunayPoly,2)
        poly = delaunayPoly[:,ind];

        locPointA = extBaseLattice[poly[1],:];
        locPointB = extBaseLattice[poly[2],:];
        locPointC = extBaseLattice[poly[3],:];
        locPointD = extBaseLattice[poly[4],:];

        midMatrix[(ind-1)*6+1,:] = (locPointA+locPointB)./2;
        midMatrix[(ind-1)*6+2,:] = (locPointA+locPointC)./2;
        midMatrix[(ind-1)*6+3,:] = (locPointA+locPointD)./2;
        midMatrix[(ind-1)*6+4,:] = (locPointB+locPointC)./2;
        midMatrix[(ind-1)*6+5,:] = (locPointB+locPointD)./2;
        midMatrix[(ind-1)*6+6,:] = (locPointC+locPointD)./2;

        thirdMatrix[(ind-1)*4+1,:] = (locPointA+locPointB+locPointC)./3;
        thirdMatrix[(ind-1)*4+2,:] = (locPointA+locPointB+locPointD)./3;
        thirdMatrix[(ind-1)*4+3,:] = (locPointB+locPointC+locPointD)./3;
        thirdMatrix[(ind-1)*4+4,:] = (locPointA+locPointC+locPointD)./3;

        fourthMatrix[ind,:] = (locPointA+locPointB+locPointC+locPointD)./4;
    end

    # Add on Delaunay Tess for Surface Atoms (2D)
    midSurfMatrix = zeros(size(surfDelaunayPoly,2)*3,3);
    thirdSurfMatrix = zeros(size(surfDelaunayPoly,2)*1,3);

    for ind = 1:size(surfDelaunayPoly,2)
        poly = surfDelaunayPoly[:,ind];

        locPointA = surfaceAtoms[poly[1],:];
        locPointB = surfaceAtoms[poly[2],:];
        locPointC = surfaceAtoms[poly[3],:];

        midSurfMatrix[(ind-1)*3+1,:] = (locPointA+locPointB)./2;
        midSurfMatrix[(ind-1)*3+2,:] = (locPointA+locPointC)./2;
        midSurfMatrix[(ind-1)*3+3,:] = (locPointB+locPointC)./2;

        thirdSurfMatrix[(ind-1)*1+1,:] = (locPointA+locPointB+locPointC)./3;
    end

    midMatrix = vcat(midMatrix,midSurfMatrix);
    thirdMatrix = vcat(thirdMatrix,thirdSurfMatrix);

    # Trim tesselation data to the simulation box.
    xkey = (midMatrix[:,1].>=0) .& (midMatrix[:,1].<xSize);
    ykey = (midMatrix[:,2].>=0) .& (midMatrix[:,2].<ySize);
    zkey = (midMatrix[:,3].>=0);
    midMatrix = midMatrix[xkey .& ykey .& zkey,:];
    midMatrix = unique(midMatrix,dims=1);

    xkey = (thirdMatrix[:,1].>=0) .& (thirdMatrix[:,1].<xSize);
    ykey = (thirdMatrix[:,2].>=0) .& (thirdMatrix[:,2].<ySize);
    zkey = (thirdMatrix[:,3].>=0);
    thirdMatrix = thirdMatrix[xkey .& ykey .& zkey,:];
    thirdMatrix = unique(thirdMatrix,dims=1);

    xkey = (fourthMatrix[:,1].>=0) .& (fourthMatrix[:,1].<xSize);
    ykey = (fourthMatrix[:,2].>=0) .& (fourthMatrix[:,2].<ySize);
    zkey = (fourthMatrix[:,3].>=0);
    fourthMatrix = fourthMatrix[xkey .& ykey .& zkey,:];
    fourthMatrix = unique(fourthMatrix,dims=1);

    # Remove Any Overlaps with Base Atoms
    extBaseLattice = BaseHaf[:,1:3];
    BLxlo = extBaseLattice[extBaseLattice[:,1].<1*Unit[1],:].+[xSize,0,0]';
    BLxhi = extBaseLattice[extBaseLattice[:,1].>xSize-(1*Unit[1]),:].-[xSize,0,0]';
    extBaseLattice = [extBaseLattice;BLxlo;BLxhi];
    BLylo = extBaseLattice[extBaseLattice[:,2].<1*Unit[2],:].+[0,ySize,0]';
    BLyhi = extBaseLattice[extBaseLattice[:,2].>ySize-(1*Unit[2]),:].-[0,ySize,0]';
    extBaseLattice = [extBaseLattice;BLylo;BLyhi];

    for ind = 1:size(extBaseLattice,1)
        locAtom = extBaseLattice[ind,:];
        keyMID = sqrt.((midMatrix[:,1].-locAtom[1]).^2+(midMatrix[:,2].-locAtom[2]).^2+(midMatrix[:,3].-locAtom[3]).^2).>MinSpacing;
        keyTHIRD = sqrt.((thirdMatrix[:,1].-locAtom[1]).^2+(thirdMatrix[:,2].-locAtom[2]).^2+(thirdMatrix[:,3].-locAtom[3]).^2).>MinSpacing;
        keyFOURTH = sqrt.((fourthMatrix[:,1].-locAtom[1]).^2+(fourthMatrix[:,2].-locAtom[2]).^2+(fourthMatrix[:,3].-locAtom[3]).^2).>MinSpacing;

        midMatrix = midMatrix[keyMID,:];
        thirdMatrix = thirdMatrix[keyTHIRD,:];
        fourthMatrix = fourthMatrix[keyFOURTH,:];
    end

    return midMatrix,thirdMatrix,fourthMatrix
end

function RecalcOxygenLattice(Box,Unit,HfLattice,OxLattice,MinSpacing)
    # Recalculate the lattice of trial oxygen points.  First uses GenerateOxygenTrialPoints to generate a lattice of points based on the atom structure in the BaseLattice (HFLatticeSites).
    # Then those points are checked against the list of existing oxygen atoms.
    xSize,ySize,zSize = Box;
	BaseOxy = pyconvert(Matrix{Any},OxLattice.positions);

    # Create ghost atoms to account for periodic boundary conditions
    extBaseOxy = BaseOxy[:,1:3];
    BLxlo = extBaseOxy[extBaseOxy[:,1].<1*Unit[1],:].+[xSize,0,0]';
    BLxhi = extBaseOxy[extBaseOxy[:,1].>xSize-(1*Unit[1]),:].-[xSize,0,0]';
    extBaseOxy = [extBaseOxy;BLxlo;BLxhi];
    BLylo = extBaseOxy[extBaseOxy[:,2].<1*Unit[2],:].+[0,ySize,0]';
    BLyhi = extBaseOxy[extBaseOxy[:,2].>ySize-(1*Unit[2]),:].-[0,ySize,0]';
    extBaseOxy = [extBaseOxy;BLylo;BLyhi];

    midMatrix,thirdMatrix,fourthMatrix = GenerateOxygenTrialPoints(Box,Unit,HfLattice,MinSpacing);

    # Scan list of trial sites to remove those within mergePointCuttoff of an existing oxygen atom.
    for ind = 1:size(extBaseOxy,1)
        locPOS = extBaseOxy[ind,1:3];
        Dist2 = sqrt.((locPOS[1].-midMatrix[:,1]).^2+(locPOS[2].-midMatrix[:,2]).^2+(locPOS[3].-midMatrix[:,3]).^2);
        Dist3 = sqrt.((locPOS[1].-thirdMatrix[:,1]).^2+(locPOS[2].-thirdMatrix[:,2]).^2+(locPOS[3].-thirdMatrix[:,3]).^2);
        Dist4 = sqrt.((locPOS[1].-fourthMatrix[:,1]).^2+(locPOS[2].-fourthMatrix[:,2]).^2+(locPOS[3].-fourthMatrix[:,3]).^2);

        keep2 = Dist2.>MinSpacing;
        keep3 = Dist3.>MinSpacing;
        keep4 = Dist4.>MinSpacing;

        midMatrix = midMatrix[keep2,:];
        thirdMatrix = thirdMatrix[keep3,:];
        fourthMatrix = fourthMatrix[keep4,:];
    end

    # Combine remaining tesselation points into a single list.
    MyTrialSites = vcat(midMatrix,thirdMatrix,fourthMatrix);
    return MyTrialSites
end

function adsorbAtom(HfAtoms,OxAtoms,TrialSites,UnitCell,BoxSizes,Temp,Press,ImpactScalingFactor=0.5)
	# Converting coordinates to the correct format
	HfPositions = pyconvert(Matrix{Any},HfAtoms.positions)

	# Weighing the possible adsorption sites with respect to the depth
	LocTrialSites = TrialSites[maximum(HfPositions[:,3]) .- TrialSites[:,3] .< 0.9*UnitCell[3],:];
	SurfWeights = Weights(exp.(-((maximum(HfPositions[:,3]) .- LocTrialSites[:,3])) .* ImpactScalingFactor));
	NumSurfSites = size(SurfWeights,1);
	LocSiteNum = sample(1:NumSurfSites,SurfWeights);

	# Choosing an adsorption site and checking if it is full
	locSite = LocTrialSites[LocSiteNum,:];
	if locSite[4] == 0
		OxAtoms.append("O");
		OxAtoms.positions[-1] = locSite[1:3];
	end

	# Time increment
	MassO2, Kb = 5.3134*10^-26, 1.380649*10^-23;
	AdsRate = ((Press*100000)./sqrt(2*pi*MassO2*Kb*Temp)*(BoxSizes[1]*BoxSizes[2]*1e-20) / 1e9)^-1;
	dt = AdsRate / (1)*log(1/(rand(1)[1]));

	return OxAtoms,dt
end

function sortNeighbours(PosNei,OxAtoms,Target,UnitCell)
	# Calculate the distances from the translating atom
	#Distances = colwise(SqEuclidean(),repeat(Target',size(PosNei,1))',PosNei');
	Distances = [norm(Target .- PosNei[i,:]) for i in 1:size(PosNei,1)];

	# Filtering between the possible destinations
	Bools = (Distances .> 0.1) .& (Distances .< UnitCell[1]);
	PosNei = hcat(PosNei[:,1][Bools],PosNei[:,2][Bools],PosNei[:,3][Bools]);
	Distances = Distances[(Distances .> 0.1) .& (Distances .< UnitCell[1])];
	PossibleNeighbors = hcat(PosNei,Distances)[sortperm(hcat(PosNei,Distances)[:, 4]), :];

	return PossibleNeighbors
end

function writeDestination(Target,Location,file="coords.dat")
	open(file, "w") do io
		println(io, string(1));
		println(io, string(Target) * " " * string(Location[1]) * " " * string(Location[2]) * " " * string(Location[3]));
	end
end

function runNEB(HfAtoms,OxAtoms,BoxSizes,TargetAtom,FinalLoc,log="log.lammps",loc="dump.neb")
	# Write NEB Input Files
	ase.io.write("Atoms.dat",HfAtoms+OxAtoms,format="lammps-data",atom_style="charge");
	writeDestination(TargetAtom,FinalLoc);
	writeNEB(HfAtoms,OxAtoms,BoxSizes,TargetAtom);

	println("NEB Started");
	command = command2;
	read(command)

	logFile = readlines(log);
    logLast = split(logFile[end])[7];
	Barrier = parse(Float64,logLast);

	locFile = readlines(loc);
    locLast = split(locFile[end])[3:5];
	Location = [parse(Float64,locLast[1]),parse(Float64,locLast[2]),parse(Float64,locLast[3])]

	println("NEB Achieved");
	return Barrier,Location
end

function logBarriers(ebs)
    #log energy barriers
    eb_log = open("eb_log.txt", "a");
    println(eb_log, ebs);
end

function timeInference(dE,Temp)
	v0,Kb_ev,h_ev,time_ns = 1,8.617333262e-5,4.135667662e-15,1e9;
	prefactor = v0 * (Kb_ev * Temp / h_ev);
	return time_ns * prefactor * exp(- dE / (Kb_ev*Temp));
end

function dumpConfiguration(HfAtoms,OxAtoms)
    # Dump atomic coordinates
	ase.io.write("HfAtoms.extxyz",HfAtoms,format="extxyz")
	ase.io.write("OxAtoms.extxyz",OxAtoms,format="extxyz")
	ase.io.write("Atoms.extxyz",HfAtoms+OxAtoms,format="extxyz")
end

function writeTime(Time,OxAtoms)
    #write time log file
    NumOxy = size(pyconvert(Any,OxAtoms.positions),1);
    Times = open("time.txt", "a")
    println(Times, NumOxy)
    println(Times, Time)
end
