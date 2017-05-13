%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  GA to find protein stable conformation states
%5%
%%%
%5% CLeandro 2014


%%
clear all
clc
close all
format long
%% Set global variables
global Ext_peptide  % Exetended peptide sequence
global amino_acid   % Dictionary with Ramachandran plots
global density      % Resolution of each Ramachandran plots

%%
% Data selection
%    You must selected here the protein

%pdb_code='2FKL'; 
%pdb_code='1PEN'; 
%pdb_code='1NOT'; 
%pdb_code='1AIE';
%pdb_code='1AJJ';
pdb_code='1FXD';
%pdb_code='1PLX';
%%
% Define discritization and ga population
%
density=0.3;   % degree of dicretization in the Ramachandran plot

constrains = 1;  % 1 statistical constrain
                 % 0 without constrains
                 
prob=1;          % 1 uses probabilites
                 % 0 number of cases
                 
PopulationSize_Data = 200 % initial population

                 
starting_population = 0; % Criteria to strat population
                         % 0 - [-pi,pi] random
                         % 1 - uses expected solution
                         % 2 - uses dihedral angles on each amino acid
                         % 3 - uses contact dihedral angles
                         
ang_part = 10;           % 2*pi/ang_part - defines the amplitude for
                         % starting values of ag variables
                        
%%
% Seletion of graphic output
draw_structure=1;  % Display and manipulate 3D molecular structure
                   % 1 - on
                   % 0 - off

plot_R = 0;     % Draw Ramachandran plots
                % 1- 2D plot
                % 2- 3D plot
                
check_prot = 0; % check dihedral angles for expected solution
                % on Ramachandran plot 
                
save_img = 1;   % Save structure Ramachandran plots in eps format
                % 1 - on
                % 0 - off

plot_Solution = 1;

%% Test Winner Filter (not good)
% Parameters for Winner filter
wiener = 0;      % Uses Wiener filter
wienerWindow=40; % Wiener filter window factor and noise
wienerNoise=2;                 

%%
% Display and manipulate 3D molecular structure
if draw_structure==1
    molviewer(strcat(pdb_code,'.pdb'));
    pause
    if save_img
        print(gcf,'-dpsc2',strcat('structure_',pdb_code,'.eps'))
    end
end

% Reads the amino acid sequence
protein = importdata(strcat(pdb_code,'_seq.tsv'),',');

% Reads the dihedral angles in the reference chain
test_angles = importdata(strcat(pdb_code,'_ang.tsv'),',');

% Computes the ramachandra plot resolution 
dim=round(2*pi/density)+1;

%%
% Extend sequence with contacts between peptide
cont_ext_petide=1;
for amino = 1:length(protein)-1
    Amino1=protein{amino};
    Amino2=protein{amino+1};
    Ext_peptide{cont_ext_petide} = strcat(Amino1);
    Ext_peptide{cont_ext_petide+1} = strcat(Amino1,'_',Amino2);
    cont_ext_petide=cont_ext_petide+2;
end
Ext_peptide{cont_ext_petide} = protein{amino+1};
%%
% Expected solution
cont_ang=1;
chains=1;
Exp_solution = Dictionary;
for amino = 1:length(test_angles)
    Exp_solution(int2str(chains)) = [Exp_solution(int2str(chains)) test_angles(amino,1)];   % phi
    Exp_solution(int2str(chains)) = [Exp_solution(int2str(chains)) test_angles(amino,2)]; % psi
    if isnan(test_angles(amino,2))
        tmp= Exp_solution(int2str(chains));
        Exp_solution(int2str(chains))=tmp(2:end-1);
        chains=chains+1;
    end
end
NstableSolutions= chains-1; % Number of stable solutions

%%
% Computes the ramachandra plot for each used amino acid
amino_acid=Dictionary; % defines a dicionary with Ramachandra plots
                      % used in the given protein
max_phi=Dictionary; % max value of dihedral angle phi
max_psi=Dictionary; % max value of dihedral angle psi

wienerh=round(2*pi/(wienerWindow*density)); % wiener window                       
for amino = 1:length(Ext_peptide)
    if not(amino_acid.containsKey(Ext_peptide{amino}))
        fprintf('Computing  %s probability densities...\n',Ext_peptide{amino});
        Amino=Ext_peptide{amino};
        % Ramachandra plot inicialization
        if mod(amino,2)==1
            M=ones(round(2*pi/density)+1,round(2*pi/density)+1);
        else
            M=ones(round(2*pi/density)+1,round(2*pi/density)+1,4);
        end
        try
            A=importdata(strcat('top500/biopython_',Amino,'_dist.tsv'),',');
            % Computing Ramachandra plot
            R=round((A+pi)/density)+1;
            num=length(A);
            for i=1:num
                if mod(amino,2)==0 && constrains
                    x=R(i,2);
                    y=R(i,3);
                    if not(isnan(x)) && not(isnan(y))
                        M(x,y,1)=M(x,y,1)+1;
                    end
                    x=R(i,2);
                    y=R(i,4);
                    if not(isnan(x)) && not(isnan(y))
                        M(x,y,2)=M(x,y,2)+1;
                    end
                    x=R(i,1);
                    y=R(i,3);
                    if not(isnan(x)) && not(isnan(y))
                        M(x,y,3)=M(x,y,3)+1;
                    end
                    x=R(i,1);
                    y=R(i,4);
                    if not(isnan(x)) && not(isnan(y))
                        M(x,y,4)=M(x,y,4)+1;
                    end
                else
                    x=R(i,1);
                    y=R(i,2);
                    if not(isnan(x)) && not(isnan(y))
                       M(x,y)=M(x,y)+1;
                    end
                end
            end
            if prob==1
                 if mod(amino,2)==0
                    M(:,:,1)=1/max(max(M(:,:,1)))*M(:,:,1);
                    M(:,:,2)=1/max(max(M(:,:,2)))*M(:,:,2);
                    M(:,:,3)=1/max(max(M(:,:,3)))*M(:,:,3);
                    M(:,:,4)=1/max(max(M(:,:,4)))*M(:,:,4);
                else
                    M=1/max(max(M))*M;
                end
            end
            if wiener==1
                if mod(amino,2)==0
                    M(:,:,1)=wiener2(M(:,:,1),[wienerh wienerh],wienerNoise);
                    M(:,:,2)=wiener2(M(:,:,2),[wienerh wienerh],wienerNoise);
                    M(:,:,3)=wiener2(M(:,:,3),[wienerh wienerh],wienerNoise);
                    M(:,:,4)=wiener2(M(:,:,4),[wienerh wienerh],wienerNoise);
                else
                    M=wiener2(M,[wienerh wienerh],wienerNoise);
                end
            end
        catch
            fprintf('Erro:  %s ...\n',Ext_peptide{amino});
        end
        amino_acid(Amino)=M; % store Ramachandra Plot
        if mod(amino,2)==1
            [v,I1]=max(max(log(M')));
            max_phi(Amino)=(I1(1)-1)*density+1/2-pi;
            [v,I2]=max(max(log(M)));
            max_psi(Amino)=(I2(1)-1)*density+1/2-pi;
        else
            M1=M(:,:,1);
            M2=M(:,:,2);
            M3=M(:,:,3);
            M4=M(:,:,4);
        end
        if plot_R==1 % Plot 2D Ramachndra Plot
            if mod(amino,2)==1
                imagesc(-log(M));
                hold on
                plot(I2(1),I1(1),'w+');
                axis image;
                title(strcat('Ramachandra plot',Amino));
                if save_img
                     print(gcf,'-dpsc2',strcat('Ramach_',Amino,'.eps'))
                end
                pause
            else
                imagesc(-log(M1));
                axis image;
                title(strcat('Ramachandra plot psi_1-phi_2:',Amino));
                if save_img
                      print(gcf,'-dpsc2',strcat('Ramach_2D_1_',Amino,'.eps'))
                end
                imagesc(-log(M2));
                title(strcat('Ramachandra plot psi_1-phi_2:',Amino));
                if save_img
                        print(gcf,'-dpsc2',strcat('Ramach_2D_2_',Amino,'.eps'))
                end
                pause 
                imagesc(-log(M3));
                title(strcat('Ramachandra plot psi_1-psi_2:',Amino));
                if save_img
                        print(gcf,'-dpsc2',strcat('Ramach_2D_3_',Amino,'.eps'))
                end
                pause
                imagesc(-log(M4));
                title(strcat('Ramachandra plot phi_1-phi_2:',Amino));
                if save_img
                        print(gcf,'-dpsc2',strcat('Ramach_2D_4_',Amino,'.eps'))
                end
            end
        end
        if plot_R==2 % Plot 3D Ramachndra Plot
            if mod(amino,2)==1
                h=surf(log(M));
                view(3);
                light;
                lighting phong;
                camlight('left');
                shading interp;
                title(strcat('Ramachandra plot',Amino));
                if save_img
                        print(gcf,'-dpsc2',strcat('Ramach_3D_',Amino,'.eps'))
                end
            else
                h=surf(log(M1));
                view(3);
                light;
                lighting phong;
                camlight('left');
                shading interp;
                title(strcat('Ramachandra plot psi_1-phi_2:',Amino));
                if save_img
                        print(gcf,'-dpsc2',strcat('Ramach_3D_1_',Amino,'.eps'))
                end
                    h=surf(log(M2));
                    view(3);
                    light;
                    lighting phong;
                    camlight('left');
                    shading interp;
                    title(strcat('Ramachandra plot psi_1-psi_2:',Amino));
                    if save_img
                        print(gcf,'-dpsc2',strcat('Ramach_3D_2_',Amino,'.eps'))
                    end
                    pause
                    h=surf(log(M3));
                    view(3);
                    light;
                    lighting phong;
                    camlight('left');
                    shading interp;
                    title(strcat('Ramachandra plot phi_1-phi_2:',Amino));
                    if save_img
                        print(gcf,'-dpsc2',strcat('Ramach_3D_3_',Amino,'.eps'))
                    end
                    pause
                    h=surf(log(M4));
                    view(3);
                    light;
                    lighting phong;
                    camlight('left');
                    shading interp;
                    title(strcat('Ramachandra plot psi_1-phsi_2:',Amino));
                    if save_img
                        print(gcf,'-dpsc2',strcat('Ramach_3D_4_',Amino,'.eps'))
                    end
                    pause
                end
        end
    end
end
%%
% Start GA optimization
%

fprintf('---------------------------------\n');
fprintf('- Starting optimal conformation  \n');
fprintf('-\n');
fprintf('- \t Protein: %s\n',pdb_code);
fprintf('- \t Population: %d\n',PopulationSize_Data);
fprintf('- \t Num. Var.: %d\n',length(Ext_peptide)-1);
fprintf('- \t Num. of stable solutions: %d\n',NstableSolutions);
minE=inf;
Exp_solution_S=Dictionary; % Normalization of native conformations
for i=1:NstableSolutions
    E=ObjFun(Exp_solution(int2str(i)));
    fprintf('- \t\t Expected free energy: %4.4f\n',E);
    if E<minE
        minE=E;
        minChain=i;
    end
    % Expected solution normalization
    %
    Exp_solution_S(int2str(i))=round((Exp_solution(int2str(i))+pi)/density)+1;
end
fprintf('---------------------------------\n');

Num_iterations=1;


% Plot dihedral angles for the expected solution
% on each correspondent Ramachandra Plot
if check_prot
    for i=2:length(protein)-1
                amino=protein{i};
                M=amino_acid(amino);
                phi = test_angles(i,1); % phi
                psi = test_angles(i,2); % psi
                rx=round((phi+pi)/density)+1;
                ry=round((psi+pi)/density)+1;
                imagesc(-log(M));
                hold on
                plot(ry,rx,'w+');
                axis image;
                title(strcat('Ramachandra plot:',amino));
                if save_img
                    print(gcf,'-dpsc2',strcat('Ramach_Test_',int2str(i),'.eps'))
                end
                hold off
                pause
    end
end

%% GA initial constants
nvars= length(Ext_peptide)-1; % number of variables

lb= -pi; % lower value 
ub=  pi; % upper value

% inital population
if starting_population==0 % lower and upper values
    PopInitRange_Data=[lb;ub];
else
    PopInitRange_Data=zeros(2,nvars);
    if starting_population==1 % uses exact solution
        for i= 1;nvars
            PopInitRange_Data(:,i)= [Exp_solution(i)-pi/ang_part;Exp_solution(i)+pi/ang_part];
        end
    elseif starting_population==2 % uses stable dihedral angles of amino acids 
            j=1;
            for i=1:2:length(Ext_peptide)
                Amino=Ext_peptide{i};
                phi=max_phi(Amino);
                psi=max_psi(Amino);
                if i>1 && i<nvars
                    PopInitRange_Data(:,j)= [phi-pi/ang_part;phi+pi/ang_part]; % initial Range
                    PopInitRange_Data(:,j+1)= [psi-pi/ang_part;psi+pi/ang_part];
                    j=j+2;
                end
                if i==1
                  PopInitRange_Data(:,j)= [psi-pi/ang_part;psi+pi/ang_part];
                  j=j+1;
                end
                if i==length(Ext_peptide)
                  PopInitRange_Data(:,j)= [phi-pi/ang_part;phi+pi/ang_part];
                end
            end
   elseif starting_population==3 % uses contact dihedral
           j=1;
           for i=2:2:length(Ext_peptide)
                Amino=Ext_peptide{i};
                phi=max_phi(Amino);
                psi=max_psi(Amino);
                PopInitRange_Data(:,j)= [phi-pi/ang_part;phi+pi/ang_part]; % initial Range
                PopInitRange_Data(:,j+1)= [psi-pi/ang_part;psi+pi/ang_part];
                j=j+2;
           end
   end
end

           

%% GA options
% Start with the default options
options = gaoptimset;

% Options setting
%% Modify options setting
options = gaoptimset(options,'TolCon', 1e-8);
options = gaoptimset(options,'PopInitRange', PopInitRange_Data);
options = gaoptimset(options,'PopulationSize', PopulationSize_Data);
options = gaoptimset(options,'FitnessScalingFcn', @fitscalingprop);

options = gaoptimset(options,'SelectionFcn', {  @selectiontournament [] });
options = gaoptimset(options,'Display', 'off');
options = gaoptimset(options,'PlotFcns', {  @gaplotbestindiv @gaplotdistance @gaplotrange @gaplotscorediversity });

for iteration=1:Num_iterations
    tic;
    [x,fval,exitflag,output,population,score] = ...
        ga(@ObjFun,nvars,[],[],[],[],lb,ub,[],[],options);
    ex_time=toc;
    hold off
    
    for i=1:length(x)
        while x(i) > pi
            x(i)=x(i)-pi;
        end
        while x(i) < -pi
            x(i)=pi+x(i);
        end
        x(i)=round((x(i)+pi)/density)+1;
    end
    local_time(iteration)= ex_time;
    fprintf('Time %4.2f\n',ex_time);
    local_solution(iteration,:) = x;
    local_value(iteration)    = fval;
    minError=inf;
    for i=1:NstableSolutions
        Error=norm(x-Exp_solution_S(int2str(i)))/length(x);
        if Error<minError
            minError=Error;
            minErrorChain=i;
        end
    end
    error_ang(iteration,:)=abs(x-Exp_solution_S(int2str(minErrorChain)));
    error(iteration)=norm(x-Exp_solution_S(int2str(minErrorChain)))/length(x);
    fprintf('Energy   %4.2f\n',fval);
    fprintf('Min. L2 error for stable conformation %d: %4.2f\n',minErrorChain,error(iteration));
    %% Expected solution normalization
    figure
    bar(Exp_solution_S(int2str(minChain)));
    title('Expected solution');
    if save_img
        print(gcf,'-dpsc2',strcat('expect_solution_',pdb_code,'.eps'))
    end
    figure
    bar(x);
    title('Best Solution')
    if save_img
        print(gcf,'-dpsc2',strcat('best_solution_',pdb_code,'.eps'))
    end
    % Plot error for each dihedral angle
    figure
    bar(error_ang(iteration,:));
    title('Variable Error')
    pause
    if save_img
             print(gcf,'-dpsc2',strcat('Error_',pdb_code,int2str(iteration),'.eps'))
    end
    % histogram with error vs number
    hist(error_ang(iteration,:));
    title('Error Histogram')
    pause
    if save_img
             print(gcf,'-dpsc2',strcat('ErrorHisto_',pdb_code,int2str(iteration),'.eps'))
    end
    %%
    % Error in each attribute
    Error_Att=Dictionary;
    num_Att=Dictionary;
    L=Exp_solution_S(int2str(minErrorChain));
    for amino=2:length(Ext_peptide)-1
        phi_r=x(amino);
        psi_r=x(amino-1);
        phi_e=L(amino);
        psi_e=L(amino-1);
        residuo=Ext_peptide{amino};
        if not(Error_Att.containsKey(Ext_peptide{amino}))
            Error_Att(residuo)=(phi_r-phi_e)^2+(psi_r-psi_e)^2;
            num_Att(residuo)=1;
        else
            Error_Att(residuo)=Error_Att(residuo)+(phi_r-phi_e)^2+(psi_r-psi_e)^2;
            num_Att(residuo)=num_Att(residuo)+1;
        end    
    end
    fprintf('------- Error by Att --------------\n');
    Residues=Error_Att.Keys();
    for i=1:length(Residues)
        residuo=Residues{i};
        fprintf('|  Att: %s  num: %d   error: %4.3f\n',residuo,num_Att(residuo),1/num_Att(residuo)*sqrt(Error_Att(residuo)));
    end
    fprintf('-----------------------------------\n');
    %%
    % Plots Ramachandra for each contact and residue
    % with real solution and the computed approximation
    %
    if plot_Solution
           for i=2:length(Ext_peptide)-1
                amino=Ext_peptide{i};
                MM=amino_acid(amino);
                if mod(i,2)==0
                    M=MM(:,:,1);
                else
                    M=MM;
                end
                px=x(i-1);py=x(i);
                rx=L(i-1);ry=L(i);
                imagesc(log(1+M));
                hold on
                plot(py,px,'w+');
                plot(ry,rx,'y+');
                axis image;
                title(strcat('Ramachandra plot:',amino));
                if save_img
                    print(gcf,'-dpsc2',strcat('Ramach_',pdb_code,amino,int2str(i),'.eps'))
                end
                pause
           end
    end
    %% 
end