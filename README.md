# FoldingWIthRotamerAndGA: Protein folding with local propensity using GA (https://melroleandro.github.io/FoldingWIthRotamerAndGA)

This work presents a set of scripts used to evaluate the usefulness of rotamer libraries to predict protein folding. The most important feature of this strategy is the possibility of modelling protein folding without explicitly treating every atom in the problem. Here I try to evaluate the importance or relevance of two  empirical energy function to predict the folding. These functions are defined using potentials extracted from Ramachandran plots (Ramakrishann 1965), studies of this plots (McCammon 1984) show that they reflect the local interactions of free energy. The protein conformation is determined from a balance between local interactions (in each amino acid) and non-local ones encoded in the protein potential energy function as statistical potentials. I explored two approximations to the free energy landscape of several proteins using the Genetic Algorithm optimization solver available in the Matlab Optimization Toolbox, trying predict a stable conformation based on the proteins first base. For that several proteins (2FKL, 1PEN, 1NOT, 1FXD) were selected from the PDB, extracting its first base, the Ramachandran plot for each of its residuum were computed and used them to approximate the protein free energy landscape. Here this empirical energy is used as a fitness function in a Genetic Algorithm Optimization procedure to predict the dihedral angles for a minimum energy conformation. For evaluate the quality of this prediction, the predict dihedral angles are compared to the protein native values. I used for similarity measurement the mean square error (MSE) between the vector of predicted dihedral angles and the value of this angles in a stable conformation.

The project is described by 2 scripts in Python, 1 MatLab script and 1 Matlab function. We used python on the data prepossessing and the MatLab to GA optimization (https://github.com/MelroLeandro/FoldingWIthRotamerAndGA.git):


 -  Extract_tests.py,is locate in the folder ./test and it is used to extract the protein 1º base and dihedral angles from a pdb file format.
 - ramachandran_Res.py, is located in the folder ./top500 and it is used to generated information used in Ramachandran plots construction. 
 - Protein_stable_multistate.m is the script used to predict protein folding using GA. This scripts reads a 1º base sequences and dihedral angles from folder .  and uses data in the folder ./top500 for density plots generation.
 -  ObjFun.m is a MatLab function used by  Protein_stable_multistate.m and the GA optimizer to compute the empirical energy surface.

The Python strips have dependences from libraries numpy, scipy,  matplotlib and BioPython.

The results and the auxiliary data are distributed by tree folders . , ./test and ./top500.

1) The folder . is the location where folding prediction graphs are saved, and where the proteins 1º base file and the file with information about the dihedral data must be located.
2) The folder ./top500 is the location for 500 proteins, in pdf format, and it location of rotamer library, for each amino acid and ech sequence of two amino acids, defined by files in text format, with respectively two and four columns. 

A protein folding prediction is generated by executing the Matlab script Protein_stable_multistate.m, here we may parametrise the solver, select the protein, select the empirical energy, configure the execution. Bellow you can see the top of this script, and test it to predict the native configuration for protein '1FXD'.   
 
%%
% Data selection
%    You must selected here the protein

%pdb_code='2FKL'; 
%pdb_code='1PEN'; 
%pdb_code='1NOT'; 
%pdb_code='1AIE';
%pdb_code='1AJJ';
pdb_code='1FXD'; % <--- my selection
%pdb_code='1PLX';
%%
% Define discritization and ga population
%
density=0.3;   % degree of dicretization in the Ramachandran plot

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
starting_population = 0; % Criteria to strat population
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
