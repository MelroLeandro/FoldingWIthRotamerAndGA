---
layout: post
title:  "Protein folding perdition using GA"
subtitle: "Using the statistical microscope"
author:  Mellean
tags:   libraries ga 
category: free-energy
visualworkflow: false
---

<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

Protein folding perdition using GA
=============

The basic of GA consists at least three operations - reproduction, crossover, and mutation. In this application domain a protein conformation is encoded as a sequence
of dihedral angles representing a chromosome. Through the application of selection operations to initial populations of chromosomes, a new generation is formed.
The initial population is here generates by the selection of random dihedral angles. The chromosomes encoding conformations with lower free energy value (our fitness
function) will be kept for breeding the next generation. The next generation then is made up of copies of conformations with high fitness which form the mating pool
for the following generation. The crossover operation mates a pair of chromosomes by randomly selecting crossover points and swapping the sequence parts. The nu-
tation randomly selects a chromosomes within the population and alters part of its dihedral angles. By applying genetic operations to an initial population of protein
conformations, a final population of lower energy conformations is formed.

The general procedure used here to predict the a stable conformation is the following:

1. Selection of a protein, in its native conformation, from the PDB.
2. Extraction of the protein 1 o base. For a 1 o base of length $n$, its conformations is encoded in a sequence of $2 ∗ n − 1$ dihedral angles used as chromosome in
the GA optimizer.
3. Generation of a rotamer library using probability densities of dihedral angles, one for each amino acid present in the protein, and for pairs of dihedral angles
in consecutive pairs of residues in the side-chain.
4. Set the GA options.
5. Execute the GA implementation.
6. Evaluate the similarity between the predicted lower energy conformations to the protein native conformation.

In each test the quality of the predicted lower energy conformation is evaluated using mean square error (MSR) between the vector of dihedral angles, for the initial 
protein conformation, and vector of dihedral angles defining the produced lower energy chromosome.

## Population

Here the chromosome is defined as a vector of dihedral angles, and its length is dependent of the selected protein. Due to restriction in the available computation
power the presented results were generated with initial population of 100 or 200
conformations. A best solution is to selected a population of n × v chromosomes,
where b is the number of bins used in the angular discritization and v is the number
of variables. For an protein defined by n residues, we have v = 2n − 1 variables or
dihedral angles. The population type is a double vector, and the initial population
is here created at random using a uniform distribution of angles in [−π, π].

## Fitness scaling

In the MatLab scaling function specifies the function that performs the scaling. In
this work the scaling function Top is used. It scales the individuals with the highest
fitness values equally.

## Reproduction

Reproduction options determine how the genetic algorithm creates children at each
new generation. Here the Elite count is set to 5 individuals that are guaranteed to
survive to the next generation. Crossover fraction specifies the fraction of the next
generation that crossover produces. Mutation produces the remaining individuals
in the next generation. The Crossover fraction was set in this work to 80%.

## Mutation

Mutation functions make small random changes in the conformations in the population, which provide genetic diversity and enable the genetic algorithm to search a
broader space. The algorithm selects a fraction of the dihedral angles for mutation,
where each angle has the same probability of being mutated. Then algorithm replaces each selected angle by a random angle selected uniformly from its domains.
Here a angle have a probability of 0.01 of being mutated.

## Crossover

Crossover combines two protein conformations, or parents, to form a new conformation, or child, for the next generation. For crossover it was selected the Scattered
method. It creates a random binary vector and selects the genes where the vector
is a 1 from the first parent, and the genes where the vector is a 0 from the second
parent, and combines the genes to form the child.

## Stopping criteria

Stopping criteria determines what causes the algorithm to terminate. The maximum
number of iterations was set to 100. If the cumulative change in the fitness function
value over Stall generations is less than $1e − 8$, the algorithm stops.

I tested the GA optimizer with different options, sets of subpopulations and
parameters, however the presented ones seems have the best performances with
short proteins.

## Example

A protein folding prediction is generated by executing the Matlab script Protein_stable_multistate.m, 
here we may parametrize the solver, select the protein, select the empirical energy, configure the execution. 
Bellow you can see the top of this script, and test it to predict the native configuration for protein ‘1FXD’.


    %%
    % Data selection
    %    You must selected here the protein


    pdb_code='1FXD'; % <--- my selection
    
    %%
    % Define discretization and ga population
    %
    density=0.3;   % degree of discretization in the Ramachandran plot

    % Selection 1 if you want dynamic constrains imposed
    % by density plots for sequences of two amino acids

    constrains = 1;  % 1 statistical constrain
                     % 0 without constrains

    % Selection of your empirical energy function
    %                 
    prob=1;          % 1 uses the probabilistic energy
                     % 0 uses the statistical energy

    % Number of random conformations for GA
    %                 
    PopulationSize_Data = 200 % GA initial population

    % Different conformation initializations  
    %                 
    starting_population = 0; % Criteria to start population
                             % 0 - [-pi,pi] random
                             % 1 - uses expected solution
                             % 2 - uses dihedral angles on each amino acid
                             % 3 - uses contact dihedral angles

    %%
    % Selection of graphic output
    draw_structure=1;  % Display and manipulate 3D molecular structure
                       % 1 - on
                       % 0 - off

    plot_R = 0;     % Draw Ramachandran plots
                    % 1- 2D plot
                    % 2- 3D plot

    check_prot = 0; % check dihedral angles for expected solution
                    % on Ramachandran plot

    save_img = 1;   % Save graphics in eps format
                    % 1 - on
                    % 0 - off

    plot_Solution = 1; % Display the Ramachandran plot with the
                       % predictions and native values.


