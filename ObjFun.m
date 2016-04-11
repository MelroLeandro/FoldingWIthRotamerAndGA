function val = ObjFun(x)
% Object function used for GA optimization
%
%
global Ext_peptide  % Exetended peptide sequence
global amino_acid   % Dictionary with Ramachandran plots
global density      % Resolution of each Ramachandran plots

m=length(Ext_peptide);

% Genome normatization
%
for i=1:length(x)
    while x(i) > pi
        x(i)=x(i)-pi;
    end
    while x(i) < -pi
        x(i)=pi+x(i);
    end
    x(i)=round((x(i)+pi)/density)+1;
end

% Free energy evaluation
%
%
val=0;
for i=2:m-1
    amino_contact=Ext_peptide{i};
    Rama=amino_acid(amino_contact); 
    if mod(i,2)==1
       % internal forces in the amino acis 
       val=val-log(Rama(x(i-1),x(i)));
    else
       % contact of type 1
       val=val-log(Rama(x(i-1),x(i),1));
       % contact of type 2
       if i<m-1
            val=val-log(Rama(x(i-1),x(i+1),2));
       end
       % contact of type 3
       if i>2
            val=val-log(Rama(x(i-2),x(i),3));
       end
       % contacy of type 4
       if i<m-2
            val=val-0.5*log(Rama(x(i-1),x(i+2),4));
       end
    end
end

end

