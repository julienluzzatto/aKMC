#!/usr/bin/env julia

# Auxiliary packages
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

# Auxiliary functions
include("LMP_EnergyEval4.jl")
include("AuxiliaryFcts.jl")

# Main function
function main()
    Run_KMC(2400,1.013,10000);
end

function Run_KMC(Temp, Press, MaxKMCmoves, ImportAtoms=false)

	# Lattice parameters (Angstroms)
    DefaultDimensions = [5,5,8]; #[6,6,20];
    HfOxygenBondDist = 2.2;
    MinOxySpacing = 1.8; #2.4;

    # Lattice generation
	if ImportAtoms == false
		println("Generating Initial Configuration")
        HfAtoms,OxAtoms = makeExistingConfiguration(DefaultDimensions);
    else
        println("Reading Initial Configuration")
        HfAtoms,OxAtoms = loadExistingConfiguration();
    end
    UnitCell,BoxSizes = getDimensions(DefaultDimensions);

    # Initial minimization
    println("Run Initial Minimization")
    HfAtoms,OxAtoms,BoxSizes = runMinimization(HfAtoms,OxAtoms,BoxSizes);
    println("Initial Configuration Obtained")

    # Begin simulation
    MoveCounter,Time = 0,0;
    ReferenceRate = timeInference(0,Temp);

    while MoveCounter < MaxKMCmoves

        # Regular minimization
        if MoveCounter % 10 == 0
            writeTime(Time,OxAtoms);
            ase.io.write("Atoms_$MoveCounter.extxyz",HfAtoms+OxAtoms,format="extxyz");

            println("----- Minimize Coordinates -----")
            HfAtoms,OxAtoms,BoxSizes = runMinimization(HfAtoms,OxAtoms,BoxSizes);
        end

        # Generate tessellation sites
        OxAtPositions = pyconvert(Matrix{Any},OxAtoms.positions);
        OxyFixedSites = hcat(OxAtPositions, 1 .* ones(size(OxAtPositions,1)));
        TessellaSites = RecalcOxygenLattice(BoxSizes,UnitCell,HfAtoms,OxAtoms,MinOxySpacing);
        OxyTrialSites = hcat(TessellaSites , 0 .* ones(size(TessellaSites, 1)));
        LocTrialSites = vcat(OxyFixedSites,OxyTrialSites);

        # Calculate rates
        NumHaf,NumOxy = size(pyconvert(Any,HfAtoms.positions),1),size(pyconvert(Any,OxAtoms.positions),1);
        Type1PerNS,Type2PerNS = 0.5,0.5 # KMC Event Rates

        println("Impact / Translation");
        println([Type1PerNS,Type2PerNS]');

        deck = Type1PerNS + Type2PerNS; # Normalize probability distribution
        draw = rand(1) * deck; # Pick which type of move

        if draw[1]<Type1PerNS # Adsorption move: an Oxygen molecule impacts surface
            println("+++++ Add up to 2 Surface Oxygen +++++")

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
            println(">>>>> Translate an Oxygen Atom <<<<<")

            # Select a random Oxygen atom
            TargetLoc = rand(1:NumOxy,1)[1];
            TargetOxy = pyconvert(Matrix{Float64},OxAtoms.positions)[TargetLoc,:];

            # Find possible destination sites
            PossibleNeighbors = TessellaSites;
            PossibleNeighbors = sortNeighbours(PossibleNeighbors,OxAtoms,TargetOxy,UnitCell);

            if ~isempty(PossibleNeighbors)
                # Restrict possible destinations
                #NumSites = size(PossibleNeighbors,1);
                NumSites = 1;
                TargetAtom = NumHaf + TargetLoc;
                EnergyBars = Array{Float64}(undef, 0, 4);

                # Calculate energy barriers
                for i=1:NumSites
                    # Determine destination coordinates
                    FinalLoc = [PossibleNeighbors[i,:][1],PossibleNeighbors[i,:][2],PossibleNeighbors[i,:][3]];

                    # Run NEB calculation
                    Barrier,Location = runNEB(HfAtoms,OxAtoms,BoxSizes,TargetAtom,FinalLoc);

                    # Add energy barrier
                    EnergyBars = vcat(EnergyBars,hcat(Barrier,FinalLoc'));
                end

                # Store energy barrier
                logBarriers(EnergyBars);
                dE = EnergyBars[1];
                k = timeInference(dE,Temp);

                OxAtoms.positions[TargetLoc,:] = Location;

                # Update time
                #if rand(1)[1] < k/ReferenceRate
                    #OxAtoms.positions[TargetLoc,:] = Location;
                #end

                dt = (- log((rand(1)[1])) / k);
                Time = Time + dt;

            else
                println("Skipping: no viable destination")
            end

        end

        dumpConfiguration(HfAtoms,OxAtoms);
        MoveCounter = MoveCounter + 1;

        println("Time / #Haf / #Oxy");
        println([Time, NumHaf, NumOxy]');
    end
end

main()
