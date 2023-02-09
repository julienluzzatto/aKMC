#!/usr/bin/env julia

#global
using ASEconvert
using LinearAlgebra
using Distributed
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

include("LMP_EnergyEval4.jl")
include("AuxiliaryFcts.jl")

#main function
function main()
    Run_KMC(2400,1.013,10000);
end

function Run_KMC(Temp, Press, MaxKMCmoves, ImportAtoms=false)

	# Lattice Parameters (Angstroms)
    DefaultDimensions = [6,6,20];
    HfOxygenBondDist = 2.2;
    MinOxySpacing = 1.65; #2.4;

    # Lattice Generation
	if ImportAtoms == false
		println("Generating Initial Configuration")
        HfAtoms,OxAtoms = makeExistingConfiguration(DefaultDimensions);
    else
        println("Reading Initial Configuration")
        HfAtoms,OxAtoms = loadExistingConfiguration();
    end
    UnitCell,BoxSizes = getDimensions(DefaultDimensions);

    # Initial Minimization
    println("Run Initial Minimization")
    HfAtoms,OxAtoms,BoxSizes = runMinimization(HfAtoms,OxAtoms,BoxSizes)
    println("Initial Configuration Obtained")

    # Begin Simulation
    MoveCounter,Time = 0,0;

    while MoveCounter < MaxKMCmoves

        # Generate Initial OxyTrialSites
        OxAtPositions = pyconvert(Matrix{Any},OxAtoms.positions);
        OxyFixedSites = hcat(OxAtPositions, 1 .* ones(size(OxAtPositions,1)));
        TessellaSites = RecalcOxygenLattice(BoxSizes,UnitCell,HfAtoms,OxAtoms,MinOxySpacing);
        OxyTrialSites = hcat(TessellaSites , 0 .* ones(size(TessellaSites, 1)));
        LocTrialSites = vcat(OxyFixedSites,OxyTrialSites);

        # Calculate Rates
        NumHaf,NumOxy = size(pyconvert(Any,HfAtoms.positions),1),size(pyconvert(Any,OxAtoms.positions),1);
        Type1PerNS,Type2PerNS = 0.5,0.5 # KMC Event Rates

        println("Impact / Translation");
        println([Type1PerNS,Type2PerNS]');

        deck = Type1PerNS + Type2PerNS; # Normalize probability distribution
        draw = rand(1) * deck; # Pick which type of move

        if draw[1]<Type1PerNS # Adsorption move: an Oxygen molecule impacts surface
            println("+++++++++++++++ Add up to 2 Surface Oxygen")

            # Add first atom to surface
            OxAtoms, dt = adsorbAtom(HfAtoms,OxAtoms,LocTrialSites,UnitCell,BoxSizes,Temp,Press);

            # Recalculate trial sites
            OxAtPositions = pyconvert(Matrix{Any},OxAtoms.positions);
            OxyFixedSites = hcat(OxAtPositions, 1 .* ones(size(OxAtPositions,1)));
            TessellaSites = RecalcOxygenLattice(BoxSizes,UnitCell,HfAtoms,OxAtoms,MinOxySpacing);
            OxyTrialSites = hcat(TessellaSites, 0 .* ones(size(TessellaSites,1)));
            LocTrialSites = vcat(OxyFixedSites,OxyTrialSites);

            # Add second atom to surface
            OxAtoms, dt = adsorbAtom(HfAtoms,OxAtoms,LocTrialSites,UnitCell,BoxSizes,Temp,Press);

            # Advance time
            Time = Time + dt;

        elseif draw[1]<Type1PerNS+Type2PerNS && NumOxy>0
            println(">>>>>>>>>>>>>>> Translate an Oxygen Atom")

            try
            # Select a Random Oxygen Atom
            TargetLoc = rand(1:NumOxy,1)[1];
            TargetOxy = OxAtoms.positions[TargetLoc];

            # Find Possible Destination Sites
            PossibleNeighbors = OxyTrialSites;
            PossibleNeighbors = sortNeighbours(PossibleNeighbors,OxAtoms,TargetLoc,UnitCell);

            if ~isempty(PossibleNeighbors)
                # Restrict Possible Destinations
                #NumSites = size(PossibleNeighbors,1);
                NumSites = 1
                TargetAtom = NumHaf + TargetLoc;
                EnergyBars = Array{Float64}(undef, 0, 4);

                # Calculate Energy Barriers
                for i=1:NumSites
                    # Determine Destination Coordinates
                    FinalLoc = [PossibleNeighbors[i,:][1],PossibleNeighbors[i,:][2],PossibleNeighbors[i,:][3]];

                    # Run NEB Calculation
                    Haf,Oxy,BoxSizes,Barrier = runNEB(HfAtoms,OxAtoms,BoxSizes,TargetAtom,FinalLoc)

                    # Add Energy Barrier
                    EnergyBars = vcat(EnergyBars,Barrier);
                end

                # Store Energy Barrier
                logBarriers(EnergyBars);
                dE = EnergyBars[1];
                k = timeInference(dE,Temp);

                # Update Time
                if dE < 5.35
                    dt = (- log((rand(1)[1])) / k);
                    Time = Time + dt;
                end

            else
                println("Skipping: no viable destination")
            end

            catch e
                showerror(stdout, e);

                translation_error = open("trans_error.txt", "a");
                println(translation_error, e);
            end
        end

        if MoveCounter % 10 == 0
            writeTime(Time,OLatticeSites);
            ase.io.write("Atoms_$MoveCounter.extxyz",HfAtoms+OxAtoms,format="extxyz")
        end
        if MoveCounter % 20 == 0
            println("----- Minimize Coordinates -----")
            HfAtoms,OxAtoms,BoxSizes = runMinimization(HfAtoms,OxAtoms,BoxSizes)
        end

        OxyTrialSites = RecalcOxygenLattice(SimDim,UnitCellSize,HFLatticeSites,OLatticeSites,MinOxySpacing);
        dumpConfiguration(HfAtoms,OxAtoms)
        MoveCounter = MoveCounter + 1;

        println("Time / #Haf / #Oxy");
        println([Time, NumHaf, NumOxy]');

        # Physical Constants (Temp in K and O partial pressure in bar)
        #Kb = 1.380649*10^-23; # Boltzmann constant [J/K]
        #MassO2 = 5.3134*10^-26; # mass of O2 [kg]
        # KMC Event Rates
        #AdsRate = ((Press*100000)./sqrt(2*pi*MassO2*Kb*Temp)*(SimDim[1]*SimDim[2]*1e-20) / 1e9)^-1; #Impact Rate estimated by molecular impingement rate (Ideal Gas)
        #TrsRate = 38.5952545516548; # Expected time for atom translation move [nanoseconds per atom]
        #ImpactScalingFactor=1.681903628276755/UnitCell[3]; # Impact strength (low numbers increase impact depth)

    end
end

main()
