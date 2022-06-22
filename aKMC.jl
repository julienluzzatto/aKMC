#!/usr/bin/env julia

using Interpolations
using UnicodePlots
using DelimitedFiles
using CMAEvolutionStrategy
using Glob
using StatsBase
using MiniQhull
using GeometricalPredicates
using LoopVectorization
using PyCall

include("LMP_EnergyEval3.jl")

## Main Function:
function main()
    KMC_param=[38.5952545516548, 1.681903628276755]; #29
    Run_KMC(2000,1.013,KMC_param,1100)
end

## Parameter Optimization Target Function
function my_MSE_Eval(KMCparams_mult)

end

function loadExistingConfiguration(BaseLatticeFile,BaseOxyFile,BaseBoxFile)
	println("Load Configuration")
	HFLatticeSites = reshape([],0,7);
    OLatticeSites = reshape([],0,7);    
	boxparams = reshape([],0,3);
	locparams=[[] [] []]
    for line in readlines(BaseLatticeFile)
        atom = reshape([parse(Float64, x) for x in split(line, "\t")],1,7)
        HFLatticeSites=vcat(HFLatticeSites,atom);
    end
    
    for line in readlines(BaseOxyFile)
        atom = reshape([parse(Float64, x) for x in split(line, "\t")],1,7)
        OLatticeSites=vcat(OLatticeSites,atom);
    end

    for line in readlines(BaseBoxFile)
        locparams = reshape([parse(Float64, x) for x in split(line, "\t")],1,3)
        boxparams=vcat(boxparams,locparams);
    end
	
	println("Configuration Loaded")
	return HFLatticeSites,OLatticeSites,boxparams
end


function generateBaseBCCLattice(RepX, RepY, RepZ, alpha::Any=3.5416)
    #alpha=3.5416;	#lattice parameters
    HFLatticeSites=reshape([],0,6);
    HFUnitCellSites=[3 0 0 0 0 1;
                 3 0 alpha/2  alpha/2  alpha/2  1];

    HFNumCellSites=size(HFUnitCellSites,1);

    for ind=0:RepX-1
        for jnd=0:RepY-1
            for knd=0:RepZ-1

                LocUnitCell=HFUnitCellSites+repeat([0 0 alpha*ind alpha*jnd alpha*knd 0],HFNumCellSites,1);
                #display(LocUnitCell);
                HFLatticeSites=vcat(HFLatticeSites,LocUnitCell);
            end
        end
    end

    HFLatticeSites=hcat(collect(1:size(HFLatticeSites,1)),HFLatticeSites);

    FROZENkey=broadcast(<,HFLatticeSites[:,6],alpha);
    NVTkey=(HFLatticeSites[:,6].<3*alpha);
    HFLatticeSites[NVTkey,2].=2;
    HFLatticeSites[FROZENkey,2].=1;

    return HFLatticeSites
end

function GenerateOxygenTrialPoints(SimDim,UnitCellSize, BaseLattice::Matrix=HFLatticeSites, mergePointCuttoff::Float64=.66)
    xSize=SimDim[1];
    ySize=SimDim[2];
    zSize=SimDim[3];

    # Create ghost atoms to account for periodic boundary conditions
    extBaseLattice=BaseLattice[:,4:6];
    BLxlo=extBaseLattice[extBaseLattice[:,1].<4*UnitCellSize[1],:].+[xSize,0,0]';
    BLxhi=extBaseLattice[extBaseLattice[:,1].>xSize-(4*UnitCellSize[1]),:].-[xSize,0,0]';
    extBaseLattice=[extBaseLattice;BLxlo;BLxhi];
    BLylo=extBaseLattice[extBaseLattice[:,2].<4*UnitCellSize[2],:].+[0,ySize,0]';
    BLyhi=extBaseLattice[extBaseLattice[:,2].>ySize-(4*UnitCellSize[2]),:].-[0,ySize,0]';
    extBaseLattice=[extBaseLattice;BLylo;BLyhi];
    BLzlo=extBaseLattice[(extBaseLattice[:,3].<3*UnitCellSize[3]) .& (extBaseLattice[:,3].>0),:].*[1,1,-1]';
    extBaseLattice=[extBaseLattice;BLzlo];

    numBasePoints=size(extBaseLattice,1);


    #Account for surface roughness
    #ID surface atoms
    function nlargest(v, n; rev=true)
        result = falses(size(v))
        result[partialsortperm(v, 1:n; rev=rev)] .= true
        return result
    end
	RepZ=round(zSize/UnitCellSize[3]);
    KEY_upperExtBL=nlargest(extBaseLattice[:,3],minimum([round(Int,numBasePoints*2/RepZ),numBasePoints]));

    Index_upperExtBL=collect(1:1:numBasePoints)


    Ind_surfaceAtoms=Index_upperExtBL[KEY_upperExtBL];
    surfaceAtoms=extBaseLattice[KEY_upperExtBL,:];
    keepAtom=zeros(Int,size(Ind_surfaceAtoms,1));

    for ind = 1:size(Ind_surfaceAtoms,1)
        locPos=surfaceAtoms[ind,:];
        searchSpacing=maximum(UnitCellSize[1:2])*3/8;
        xdist=surfaceAtoms[:,1].-locPos[1];
        ydist=surfaceAtoms[:,2].-locPos[2];
        neighKEY=(abs.(xdist).<searchSpacing).&(abs.(ydist).<searchSpacing);
        neighZheights=surfaceAtoms[neighKEY,3];

        isheighest=locPos[3]>=maximum(neighZheights); #loc is highest? true false
        keepAtom[ind]=isheighest;
    end

    surfaceAtoms=surfaceAtoms[convert.(Bool,keepAtom),:];
    Ind_surfaceAtoms=Ind_surfaceAtoms[convert.(Bool,keepAtom)];
    surfDelaunayPoly=delaunay(2,size(surfaceAtoms,1),vec(surfaceAtoms[:,1:2]'))

    #Calcluate Regular Delaunay Tess
    delaunayPoly=delaunay(3,numBasePoints,vec(extBaseLattice[:,1:3]'));
    midMatrix=zeros(size(delaunayPoly,2)*6,3);
    thirdMatrix=zeros(size(delaunayPoly,2)*4,3);
    fourthMatrix=zeros(size(delaunayPoly,2)*1,3);

    #Calculate 2,3, and 4 fold midpoints
    for ind=1:size(delaunayPoly,2)
        poly=delaunayPoly[:,ind];

        locPointA=extBaseLattice[poly[1],:];
        locPointB=extBaseLattice[poly[2],:];
        locPointC=extBaseLattice[poly[3],:];
        locPointD=extBaseLattice[poly[4],:];


        midMatrix[(ind-1)*6+1,:]=(locPointA+locPointB)./2;
        midMatrix[(ind-1)*6+2,:]=(locPointA+locPointC)./2;
        midMatrix[(ind-1)*6+3,:]=(locPointA+locPointD)./2;
        midMatrix[(ind-1)*6+4,:]=(locPointB+locPointC)./2;
        midMatrix[(ind-1)*6+5,:]=(locPointB+locPointD)./2;
        midMatrix[(ind-1)*6+6,:]=(locPointC+locPointD)./2;

        thirdMatrix[(ind-1)*4+1,:]=(locPointA+locPointB+locPointC)./3;
        thirdMatrix[(ind-1)*4+2,:]=(locPointA+locPointB+locPointD)./3;
        thirdMatrix[(ind-1)*4+3,:]=(locPointB+locPointC+locPointD)./3;
        thirdMatrix[(ind-1)*4+4,:]=(locPointA+locPointC+locPointD)./3;

        fourthMatrix[ind,:]=(locPointA+locPointB+locPointC+locPointD)./4;
    end

    # Add on Delaunay Tess for Surface Atoms (2D)
    midSurfMatrix=zeros(size(surfDelaunayPoly,2)*3,3);
    thirdSurfMatrix=zeros(size(surfDelaunayPoly,2)*1,3);

    for ind=1:size(surfDelaunayPoly,2)
        poly=surfDelaunayPoly[:,ind];

        locPointA=surfaceAtoms[poly[1],:];
        locPointB=surfaceAtoms[poly[2],:];
        locPointC=surfaceAtoms[poly[3],:];

        midSurfMatrix[(ind-1)*3+1,:]=(locPointA+locPointB)./2;
        midSurfMatrix[(ind-1)*3+2,:]=(locPointA+locPointC)./2;
        midSurfMatrix[(ind-1)*3+3,:]=(locPointB+locPointC)./2;

        thirdSurfMatrix[(ind-1)*1+1,:]=(locPointA+locPointB+locPointC)./3;

    end


    midMatrix = vcat(midMatrix,midSurfMatrix);
    thirdMatrix = vcat(thirdMatrix,thirdSurfMatrix);


    # Trim tesselation data to the simulation box.
    xkey=(midMatrix[:,1].>=0) .& (midMatrix[:,1].<xSize);
    ykey=(midMatrix[:,2].>=0) .& (midMatrix[:,2].<ySize);
    zkey=(midMatrix[:,3].>=0);
    midMatrix=midMatrix[xkey .& ykey .& zkey,:];
    midMatrix=unique(midMatrix,dims=1);

    xkey=(thirdMatrix[:,1].>=0) .& (thirdMatrix[:,1].<xSize);
    ykey=(thirdMatrix[:,2].>=0) .& (thirdMatrix[:,2].<ySize);
    zkey=(thirdMatrix[:,3].>=0);
    thirdMatrix=thirdMatrix[xkey .& ykey .& zkey,:];
    thirdMatrix=unique(thirdMatrix,dims=1);

    xkey=(fourthMatrix[:,1].>=0) .& (fourthMatrix[:,1].<xSize);
    ykey=(fourthMatrix[:,2].>=0) .& (fourthMatrix[:,2].<ySize);
    zkey=(fourthMatrix[:,3].>=0);
    fourthMatrix=fourthMatrix[xkey .& ykey .& zkey,:];
    fourthMatrix=unique(fourthMatrix,dims=1);


    # Remove Any Overlaps with Base Atoms
    extBaseLattice=BaseLattice[:,4:6];
    BLxlo=extBaseLattice[extBaseLattice[:,1].<1*UnitCellSize[1],:].+[xSize,0,0]';
    BLxhi=extBaseLattice[extBaseLattice[:,1].>xSize-(1*UnitCellSize[1]),:].-[xSize,0,0]';
    extBaseLattice=[extBaseLattice;BLxlo;BLxhi];
    BLylo=extBaseLattice[extBaseLattice[:,2].<1*UnitCellSize[2],:].+[0,ySize,0]';
    BLyhi=extBaseLattice[extBaseLattice[:,2].>ySize-(1*UnitCellSize[2]),:].-[0,ySize,0]';
    extBaseLattice=[extBaseLattice;BLylo;BLyhi];


    for ind=1:size(extBaseLattice,1)
        locAtom=extBaseLattice[ind,:];
        keyMID=sqrt.((midMatrix[:,1].-locAtom[1]).^2+(midMatrix[:,2].-locAtom[2]).^2+(midMatrix[:,3].-locAtom[3]).^2).>mergePointCuttoff;
        keyTHIRD=sqrt.((thirdMatrix[:,1].-locAtom[1]).^2+(thirdMatrix[:,2].-locAtom[2]).^2+(thirdMatrix[:,3].-locAtom[3]).^2).>mergePointCuttoff;
        keyFOURTH=sqrt.((fourthMatrix[:,1].-locAtom[1]).^2+(fourthMatrix[:,2].-locAtom[2]).^2+(fourthMatrix[:,3].-locAtom[3]).^2).>mergePointCuttoff;

        midMatrix=midMatrix[keyMID,:];
        thirdMatrix=thirdMatrix[keyTHIRD,:];
        fourthMatrix=fourthMatrix[keyFOURTH,:];
    end


    return midMatrix,thirdMatrix,fourthMatrix

end

function RecalcOxygenLattice(SimDim,UnitCellSize, BaseLattice::Matrix=HFLatticeSites, BaseOxy::Matrix=OLatticeSites, mergePointCuttoff::Float64=1.65)
# Recalculate the lattice of trial oxygen points.  First uses GenerateOxygenTrialPoints to generate a lattice of points based on the atom structure in the BaseLattice (HFLatticeSites).
# Then those points are checked against the list of existing oxygen atoms.

    xSize=SimDim[1];
    ySize=SimDim[2];
    zSize=SimDim[3];

    # Create ghost atoms to account for periodic boundary conditions
    extBaseOxy=BaseOxy[:,4:6];
    BLxlo=extBaseOxy[extBaseOxy[:,1].<1*UnitCellSize[1],:].+[xSize,0,0]';
    BLxhi=extBaseOxy[extBaseOxy[:,1].>xSize-(1*UnitCellSize[1]),:].-[xSize,0,0]';
    extBaseOxy=[extBaseOxy;BLxlo;BLxhi];
    BLylo=extBaseOxy[extBaseOxy[:,2].<1*UnitCellSize[2],:].+[0,ySize,0]';
    BLyhi=extBaseOxy[extBaseOxy[:,2].>ySize-(1*UnitCellSize[2]),:].-[0,ySize,0]';
    extBaseOxy=[extBaseOxy;BLylo;BLyhi];

    midMatrix,thirdMatrix,fourthMatrix =GenerateOxygenTrialPoints(SimDim,UnitCellSize,BaseLattice);

    # Scan list of trial sites to remove those within mergePointCuttoff of an existing oxygen atom.
    for ind=1:size(extBaseOxy,1)
        locPOS=extBaseOxy[ind,1:3];
        Dist2=sqrt.((locPOS[1].-midMatrix[:,1]).^2+(locPOS[2].-midMatrix[:,2]).^2+(locPOS[3].-midMatrix[:,3]).^2);
        Dist3=sqrt.((locPOS[1].-thirdMatrix[:,1]).^2+(locPOS[2].-thirdMatrix[:,2]).^2+(locPOS[3].-thirdMatrix[:,3]).^2);
        Dist4=sqrt.((locPOS[1].-fourthMatrix[:,1]).^2+(locPOS[2].-fourthMatrix[:,2]).^2+(locPOS[3].-fourthMatrix[:,3]).^2);

        keep2=Dist2.>mergePointCuttoff;
        keep3=Dist3.>mergePointCuttoff;
        keep4=Dist4.>mergePointCuttoff;

        midMatrix=midMatrix[keep2,:];
        thirdMatrix=thirdMatrix[keep3,:];
        fourthMatrix=fourthMatrix[keep4,:];
    end


    # Combine remaining tesselation points into a single list.
    MyTrialSites=vcat(midMatrix,thirdMatrix,fourthMatrix);
    return MyTrialSites
end

function dump_CONFIG()

    global OxyTrialSites
    global HFLatticeSites
    global OLatticeSites
    #Dump points to text files

    open("trialOxygen.dat", "w") do io
    writedlm(io, OxyTrialSites)
    end;

    open("base.dat", "w") do io
    writedlm(io, HFLatticeSites)
    end;

    open("baseOxygen.dat", "w") do io
       writedlm(io, OLatticeSites)
    end;
end

function write_traj()

    global HFLatticeSites
    global OLatticeSites
    global Time
    #Dump trajectory to text files

    nb_Hf,nb_O = size(HFLatticeSites)[1],size(OLatticeSites)[1]

    traj = open("trajectory.txt", "a")

    println(traj, nb_Hf+nb_O)
    println(traj, Time)

    count_Hf,count_O = 1,1
    while count_Hf < nb_Hf
        println(traj, join(repr.(HFLatticeSites[count_Hf,4:6]),", "))
        count_Hf = count_Hf + 1
    end;

    while count_O < nb_O
        println(traj, join(repr.(OLatticeSites[count_O,4:6]),", "))
        count_O = count_O + 1
    end;
end

function write_time()

    global Time
    global OLatticeSites

    NumbOxy=size(OLatticeSites,1);

    traj = open("time.txt", "a")

    println(traj, NumbOxy)
    println(traj, Time)

end

function Run_KMC(Temp, Press, KMCparams, MaxKMCtime)
    ImportAtoms=true;

	# Lattice Generation Parameters
    HfOxygenBondDist=1.77;   #angstroms (Covalent radius or S-orbital radius) is also approx sigma LJ parameter
    MinOxySpacing=1.65; #angstroms 	spacing used when grid searching for floating lattice points

    #Initialize HF Lattice
    global HFLatticeSites = reshape([],0,7);
    global OLatticeSites = reshape([],0,7);    
	
	if ImportAtoms==true
		println("Reading Initial Configuration")
		HFLatticeSites,OLatticeSites,boxparams=loadExistingConfiguration("baseHf_hcp.dat","baseO_hcp.dat","baseDIM_hcp.dat");
		global SimDim=boxparams[1,:];
		global UnitCellSize=boxparams[2,:];
	else
		println("Generating Initial Configuration")
		alpha=3.5416;	#lattice parameters (assuming cubic)

		RepX=6;		#lattice unit cells in X dim
		RepY=6;		#lattice unit cells in Y dim
		RepZ=10;	#lattice unit cells in Z dim

		XwallHi=alpha*RepX;
		YwallHi=alpha*RepY;
		ZwallHi=alpha*RepZ+50; #Hight of material + 50 angstrom buffer
		
		global SimDim=[XwallHi, YwallHi, ZwallHi];	#Size of simulation box in angstroms (leave at least 20 angstroms of empty space above surface).
		global UnitCellSize=[alpha, alpha, alpha]#[3.2,2.7,5.08]; #Approx unit cell size in each dimension (rectangular simulation box only).
		HFLatticeSites=generateBaseBCCLattice(RepX,RepY,RepZ);
	end

    #Set KMC parameters
		#Temp in Kelvin
		#Press = oxygen partial pressure in bar
    dataEvery=2; #Output data every # of attempted moves

    Kb=1.380649*10^-23; #Boltzmann constant [J/K]
    Kb_ev=8.617333262145*10^-5; #Boltzmann constant [J/K]

    MassO2=5.3134*10^-26; # Mass of O2 [kg]

    #KMC Event rates
    AdsRate=((Press*100000)./sqrt(2*pi*MassO2*Kb*Temp)*(SimDim[1]*SimDim[2]*1e-20) / 1e9)^-1; #Impact Rate estimated by molecular impingement rate (Ideal Gas)
    TrsRate=KMCparams[1]; # Expected time for atom translation move [nanoseconds per atom]
    ImpactScalingFactor=KMCparams[2]/UnitCellSize[3]; # Impact strength (low numbers increase impact depth)

    #MD parameters
    global MD_timestep=.0005;
    SetMD_Sims=8;
    MaxConcurrentSims=4;
    smallMinSteps=10;
    smallMDsteps=10;
    mediumMinSteps=50;
    mediumMDsteps=250;
    largeMDsteps=25000;
	surfaceDepthInCells=0.9
    #surfaceDepthInCells=1.5;
    global LMPvect=startN_LAMMPS_instances(SetMD_Sims);

	
   println("Run Initial Minimize")
   OLatticeSites,HFLatticeSites = MinimizeCoords(LMPvect,OLatticeSites,HFLatticeSites,Temp,MD_timestep,SimDim,smallMinSteps,smallMDsteps);
   println("Initial Configuration Obtained")
    #Generate Initial OxyTrialSites
    global OxyTrialSites=RecalcOxygenLattice(SimDim,UnitCellSize,HFLatticeSites,OLatticeSites,MinOxySpacing);


    ## Begin Simulation
    global MoveCounter=0; #Number of moves taken
    global Time=0; #nanoseconds
    global LastTime=0; #nanoseconds
    global OxAdsorbed=[0 0]; #Number of adsorbed oxygen atoms [time,#atoms]
    global PossibleNeighbors=[];


    global OxyTrialSites=RecalcOxygenLattice(SimDim,UnitCellSize,HFLatticeSites,OLatticeSites,MinOxySpacing);
    
    #

    while MoveCounter<1000000 #Time<MaxKMCtime
        global OLatticeSites
        global HFLatticeSites
        global OxyTrialSites
        global MoveCounter
        global Time
        global LastTime
        global OxAdsorbed
        global LMPvect
        global PossibleNeighbors
        global MD_timestep
        global SimDim

        NumbHaf=size(HFLatticeSites,1)
        NumbOxy=size(OLatticeSites,1);
        Type1PerNS=1/AdsRate;			#Current O2 impact rate
        Type2PerNS=1/(TrsRate/NumbOxy);	#Current Oxygen translation move rate
        Type3PerNS=1/AdsRate/400;		#Current Probability of Running a short MD segment

        display([Type1PerNS,Type2PerNS,Type3PerNS]')

        FPdeck=Type1PerNS+Type2PerNS+Type3PerNS;	#Normalize probability distribution
        draw=rand(1)*FPdeck;			#Pick which type of move
        #display(draw)

        if draw[1]<Type1PerNS			#Oxygen molecule impacts surface
            #Add atom to surface
            println("Add up to 2 Surface Oxygen")
            #Advance time
            dt = AdsRate/(1)*log(1/(rand(1)[1]))
            Time=Time+dt;

            LocTrialSites=hcat(OxyTrialSites,zeros(size(OxyTrialSites,1),1));
            LocTrialSites=vcat(OLatticeSites[:,4:7],LocTrialSites);
            LocTrialSites=LocTrialSites[maximum(HFLatticeSites[:,6]).-LocTrialSites[:,3].<0.9*UnitCellSize[3],:];

            SurfWeights=Weights( exp.(-((maximum(HFLatticeSites[:,6]).-LocTrialSites[:,3])).*ImpactScalingFactor) );
            NumSurfSites=size(SurfWeights,1);
            LocSiteNum=sample(1:NumSurfSites,SurfWeights,2);

            locSite1=LocTrialSites[LocSiteNum[1],:];
            locSite2=LocTrialSites[LocSiteNum[2],:];

            if locSite1[4] == 0
                OLatticeSites=vcat(OLatticeSites,[size(OLatticeSites,1)+1 4 0 locSite1[1] locSite1[2] locSite1[3] 1]);
            end
            if locSite2[4] == 0
                oxysep=sqrt((locSite1[1]-locSite2[1])^2+(locSite1[2]-locSite2[2])^2+(locSite1[3]-locSite2[3])^2);
                if oxysep > MinOxySpacing
                    OLatticeSites=vcat(OLatticeSites,[size(OLatticeSites,1)+1 4 0 locSite2[1] locSite2[2] locSite2[3] 1]);
                end
            end
            if (locSite1[4] == 0) || (locSite2[4] == 0)
                OLatticeSites,HFLatticeSites = MinimizeCoords(LMPvect,OLatticeSites,HFLatticeSites,Temp,MD_timestep,SimDim,smallMinSteps,smallMDsteps);
                OxyTrialSites=RecalcOxygenLattice(SimDim,UnitCellSize,HFLatticeSites,OLatticeSites,MinOxySpacing);
            end


        elseif draw[1]<Type1PerNS+Type2PerNS && NumbOxy>0
            println("Translate an Oxygen Atom")

            try

            indi = rand(1:size(OLatticeSites,1),1)
            LocOxy=OLatticeSites[indi,:];	#Select a random oxygen atom

            ## vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Begin Translation Move

            PossibleNeighbors=OxyTrialSites;

            # Check X Periodic Boundary Conditions
            Xlo=LocOxy[4]-UnitCellSize[1]+1.0;
            Xhi=LocOxy[4]+UnitCellSize[1]-1.0;
            if Xlo<0
                Xlo2=Xlo+SimDim[1];
                PossibleNeighbors=PossibleNeighbors[(PossibleNeighbors[:,1].>Xlo2) .| (PossibleNeighbors[:,1].<Xhi),:];
            elseif Xhi>SimDim[1]
                Xhi2=Xhi-SimDim[1];
                PossibleNeighbors=PossibleNeighbors[(PossibleNeighbors[:,1].>Xlo) .| (PossibleNeighbors[:,1].<Xhi2),:];
            else
                PossibleNeighbors=PossibleNeighbors[(PossibleNeighbors[:,1].<Xhi) .& (PossibleNeighbors[:,1].>Xlo),:];
            end

            # Check Y Periodic Boundary Conditions
            Ylo=LocOxy[5]-UnitCellSize[2]+1.0;
            Yhi=LocOxy[5]+UnitCellSize[2]-1.0;
            if Ylo<0

                Ylo2=Ylo+SimDim[2];
                PossibleNeighbors=PossibleNeighbors[(PossibleNeighbors[:,2].>Ylo2) .| (PossibleNeighbors[:,2].<Yhi),:];
            elseif Yhi>SimDim[2]

                Yhi2=Yhi-SimDim[2];
                PossibleNeighbors=PossibleNeighbors[(PossibleNeighbors[:,2].>Ylo) .| (PossibleNeighbors[:,2].<Yhi2),:];
            else

                PossibleNeighbors=PossibleNeighbors[(PossibleNeighbors[:,2].<Yhi) .& (PossibleNeighbors[:,2].>Ylo),:];
            end

            # Z dim (no PBC)
            Zlo=LocOxy[6]-UnitCellSize[3]+1.0;
            Zhi=LocOxy[6]+UnitCellSize[3]-1.0;
            PossibleNeighbors=PossibleNeighbors[(PossibleNeighbors[:,3].<Zhi) .& (PossibleNeighbors[:,3].>Zlo),:];

            if ~isempty(PossibleNeighbors)
                #start with original site
                PossibleNeighbors=hcat(PossibleNeighbors,zeros(size(PossibleNeighbors,1),1));
                # sub-sample possible neighbor list
                if size(PossibleNeighbors,1)>MaxConcurrentSims-1
                    PossibleNeighbors=PossibleNeighbors[sample(1:size(PossibleNeighbors,1),MaxConcurrentSims-1, replace=false),:];
                end
                PossibleNeighbors=vcat(LocOxy[4:7]',PossibleNeighbors);
                #display(PossibleNeighbors)

                NumSites=size(PossibleNeighbors,1);
				#display(LocOxy)
                #OLatticeSites=OLatticeSites[LocOxy[1].!=OLatticeSites[:,1],:]; #Remove Moving Atom From List

                ##################################################################

                dump_CONFIG()

                open("posNeighbors.dat", "w") do io
                writedlm(io, PossibleNeighbors)
                end;

                py"""

                def dimer_search(oxygen_tag,xdim,ydim,zdim,dimer_searches=5):

                    #auxiliary packages
                    import os
                    import math
                    import copy
                    import time
                    import numpy as np
                    import scipy as sp

                    from ase import Atoms, Atom
                    from ase import io

                    from ase.io.trajectory import Trajectory
                    from ase.visualize import view
                    from ase.build import fcc100, add_adsorbate
                    from ase.constraints import FixAtoms

                    from ase.optimize import BFGS
                    from ase.dimer import DimerControl, MinModeAtoms, MinModeTranslate
                    from ase.calculators.lammpsrun import LAMMPS
                    from ase.calculators.lammpslib import LAMMPSlib

                    #environment variables

                    os.environ['LAMMPS_PATH'] = "/export/apps/lammps-29Sep21/"
                    os.environ['LAMMPS_POTENTIALS'] = "/export/apps/lammps-29Sep21/potentials/"
                    os.environ['LAMMPS_COMMAND'] = "/export/apps/lammps-29Sep21/build/lmp"
                    os.environ['WORKDIR'] = "/home/jluzz/aKMC-v3_2"
                    os.environ['ASE_LAMMPSRUN_COMMAND'] = "/export/apps/lammps-29Sep21/build/lmp"

                    #data post-processing
                    Hffile = "base.dat"
                    Ofile = "baseOxygen.dat"
                    #Otrial = "trialOxygen.dat"
                    Otrial = "posNeighbors.dat"

                    lines_Hf = []
                    with open(Hffile) as f:
                        for line in f.readlines():
                            line = line.split()
                            lines_Hf.append(['Hf',float(line[3]),float(line[4]),float(line[5])])
                    f = open("baseHf_nice.xyz", "w")
                    f.write(str(len(lines_Hf)))
                    f.write('\n')
                    f.write("")
                    f.write('\n')
                    for line in lines_Hf:
                        for el in line:
                            f.write(str(el))
                            f.write(" ")
                        f.write('\n')
                    f.close()
                    Hf = io.read('baseHf_nice.xyz')

                    lines_O = []
                    with open(Ofile) as f:
                        for line in f.readlines():
                            line = line.split()
                            lines_O.append(['O',float(line[3]),float(line[4]),float(line[5])])
                    f = open("baseOxygen_nice.xyz", "w")
                    f.write(str(len(lines_O)))
                    f.write('\n')
                    f.write("")
                    f.write('\n')
                    for line in lines_O:
                        for el in line:
                            f.write(str(el))
                            f.write(" ")
                        f.write('\n')
                    f.close()
                    Oxy = io.read('baseOxygen_nice.xyz')

                    lines_O_tr = []
                    with open(Otrial) as f:
                        for line in f.readlines():
                            line = line.split()
                            lines_O_tr.append(['O',float(line[0]),float(line[1]),float(line[2])])
                    f = open("baseOxygenTr_nice.xyz", "w")
                    f.write(str(len(lines_O_tr)))
                    f.write('\n')
                    f.write("")
                    f.write('\n')
                    for line in lines_O_tr:
                        for el in line:
                            f.write(str(el))
                            f.write(" ")
                        f.write('\n')
                    f.close()
                    OxyTr = io.read('baseOxygenTr_nice.xyz')

                    #auxiliary dimer method functions

                    #boltzmann acceptance criterion
                    def kmc_boltz(dE,v0=1,Kb_ev=8.617333262145e-5,Temp=2400,upper_bound=1):
                        return v0*np.exp(-dE/(Kb_ev*Temp),dtype=np.float128) / upper_bound

                    #construction of the structure
                    both = Hf + Oxy

                    both.set_cell([xdim,ydim,zdim])
                    both.set_pbc([True, True, False])
                    target = int(len(Hf) + oxygen_tag - 1)
                    #print(target)

                    n = len(both)
                    Kb = 1.380649e-23
                    p0 = time.time()

                    #initial position
                    r0 = both.positions[target]
                    r0b = copy.deepcopy(r0)
                    #print(r0)

                    #setting the mask
                    mask = [True] * len(both)
                    mask[target] = False
                    constraint = FixAtoms(mask=mask)
                    both.set_constraint(constraint)
                    #print(mask)

                    #auxiliary parameters
                    parameters = {'pair_style': 'comb',
                                  'pair_coeff': ['ffield.comb O Hf']}
                    #lammps = LAMMPS(parameters=parameters, files=files)
                    lammps = LAMMPS()

                    #setting the calculator
                    both.calc = lammps
                    #both.calc = LennardJones(sigma=sig, epsilon=eps)
                    e0 = both.get_potential_energy()
                    #print('e0,r0: ({},{})'.format(str(e0),str(r0)))

                    # dimer search method
                    def dimer_search():

                        #trajectory log
                        #traj = Trajectory('dimer_both.traj', 'w', both)
                        #traj.write()
                        d_mask = [not i for i in mask]

                        #setting the dimer up
                        d_control = DimerControl(#initial_eigenmode_method='displacement',
                                                 #displacement_method='vector',
                                                 initial_eigenmode_method='gauss',
                                                 displacement_method='gauss',
                                                 logfile=None,
                                                 mask=d_mask)
                        d_atoms = MinModeAtoms(both, d_control)

                        #displace the dimer
                        #displacement_vector = np.zeros((n, 3))
                        #displacement_vector[target] = [0,0,-0.0001]
                        gauss_std = 0.005
                        fmax = 0.05
                        d_atoms.displace(gauss_std=gauss_std) #displacement_vector=displacement_vector
                        dim_rlx = MinModeTranslate(d_atoms,
                                               #trajectory=traj,
                                               #logfile='logfile.txt'
                                               logfile=None)
                        dim_rlx.run(fmax=fmax)
                        #open('logfile.txt', 'w').close()

                        #trajectory
                        trajectory = False
                        if trajectory:
                            loadtraj = Trajectory('dimer_both.traj')
                            logtraj = []
                            for atoms in loadtraj:
                                logtraj.append(atoms)

                        #final position
                        #print(d_atoms.get_positions()[target])

                        #evaluate the energy barrier
                        #eb = both.get_potential_energy()
                        #diff = eb - e0
                        rb = copy.deepcopy(both.positions[target])

                        return rb

                    #perform a finite amount of dimer searches
                    saddle_points = []
                    #search_time = time.time()
                    for _ in range(dimer_searches):
                        #saddle_trial,saddle_trial_round = dimer_search()
                        saddle_trial = dimer_search()
                        #if not any((saddle_trial_round == x).all() for x in saddle_points_round):
                        saddle_points.append(saddle_trial)
                        both.positions[target] = r0b
                    m = len(saddle_points)
                    #print(m)

                    #pick a saddle point
                    j = np.random.randint(len(saddle_points))
                    rb = saddle_points[j]
                    both.positions[target] = rb
                    eb = both.get_potential_energy()
                    diff = eb - e0
                    #print(diff)

                    #iterate towards next position
                    #cutoff = 3
                    #ef_tmp,rf_tmp = np.inf,None
                    #for point in OxyTr.get_positions():
                        #if sp.spatial.distance.euclidean(r0b,point) < cutoff:
                        #both.positions[target] = point
                        #if both.get_potential_energy() < ef_tmp:
                            #rf_tmp = point
                        #both.positions[target] = rb

                    displacement = rb - r0b
                    rf_tmp = r0b + 1.025 * displacement

                    #accept or reject translation
                    accept = False
                    both.positions[target] = r0b
                    print(kmc_boltz(diff))
                    if np.random.random() < kmc_boltz(diff):
                        both.positions[target] = rf_tmp
                        rf = rf_tmp
                        accept = True
                    else:
                        rf = r0
                        both.positions[target] = r0b
                    #print('ef,rf,fo: ({},{})'.format(str(ef),str(rf),str(fo)))
                    print('Accepted: ' + str(accept) + '<<<<<<<<<<<<<<<<')

                    #minimize single atom
                    #dyn = BFGS(both)
                    #dyn.run(fmax=0.05)

                    return diff,rf
                """

                diff,rf = py"dimer_search"(indi,SimDim[1],SimDim[2],SimDim[3],1)

                #OLatticeSites=vcat(OLatticeSites,[LocOxy[1] 4 0 PossibleNeighbors[moveAccepted,1:3]' 1]);
                OLatticeSites=OLatticeSites[LocOxy[1].!=OLatticeSites[:,1],:]; #Remove Moving Atom From List
                r0bx, r0by, r0bz = rf
                OLatticeSites=vcat(OLatticeSites,[LocOxy[1] 4 0 r0bx r0by r0bz 1]);
                dt = 1/(TrsRate)*log(1/(rand(1)[1]))
                Time=Time+dt;

                ##################################################################

                if MoveCounter % 8 == 0
                    println("Minimize Coords<<<<<<<<<<<<<<<")
                    OLatticeSites,HFLatticeSites = MinimizeCoords(LMPvect,OLatticeSites,HFLatticeSites,Temp,MD_timestep,SimDim,smallMinSteps,smallMDsteps); # Minimize Coordinates
                end
                #println(Time)

                OxyTrialSites=RecalcOxygenLattice(SimDim,UnitCellSize,HFLatticeSites,OLatticeSites,MinOxySpacing);

            else
                println("skipping: no viable destinations")
            end

            catch e
                showerror(stdout, e)

				println("<><><><><><> Attempt Recovery From Failed Dimer Search <><><><><><>")

            end

            ## ^^^^^^^^^^^^^^^^^^^ End Translation Move
        else
            println("Run Short MD Simulation")
            #@time OLatticeSites,HFLatticeSites = MinimizeCoords(LMPvect,OLatticeSites,HFLatticeSites,Temp,MD_timestep,SimDim,mediumMinSteps,largeMDsteps);
            OxyTrialSites=RecalcOxygenLattice(SimDim,UnitCellSize,HFLatticeSites,OLatticeSites,MinOxySpacing);
			Time=Time+(MD_timestep*largeMDsteps*.001);
        end

        if MoveCounter % 20 == 0
            write_time()
            if isfile("log_trial.txt")
                rm("log_trial.txt")
                open("log_trial.txt","a")
            end
            if isfile("log.lammps")
                rm("log.lammps")
                open("log.lammps","a")
            end

            x_hf,y_hf,z_hf,b_hf = HFLatticeSites[:,4],HFLatticeSites[:,5],HFLatticeSites[:,6],HFLatticeSites[:,7];
			len_hf = size(x_hf)[1]
			HF_table = Array{Union{Float64,Int64,String}}(undef, len_hf, 4)
			for i in 1:len_hf
				HF_table[i,:] = ["Hf" x_hf[i] y_hf[i] z_hf[i]]
			end

            x_o,y_o,z_o,b_o = OLatticeSites[:,4],OLatticeSites[:,5],OLatticeSites[:,6],OLatticeSites[:,7];
            len_o = size(x_o)[1] + 2
            O_table = [0 " " " " " "; "#" " " " " " "]
            for i in 1:len_o-2
                j = i+2
                if b_o[i]==1
                    O_table = vcat(O_table, ["O" x_o[i] y_o[i] z_o[i]])
                end
            end
            len_O = size(O_table)[1]

            full_table = vcat(O_table,HF_table)
            full_table[1] = Int64(len_O+len_hf-2)
            writedlm("plot_xyz.$MoveCounter.xyz",full_table);
        end;

        MoveCounter=MoveCounter+1
        dump_CONFIG();
        println("Time #Haf #Oxy");
        println([Time,NumbHaf,NumbOxy]');
    end
end


main()
