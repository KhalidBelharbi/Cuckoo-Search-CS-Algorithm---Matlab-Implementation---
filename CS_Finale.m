% /**
%* Copyright - This code contains a watermark for the author
%* @author Khalid Belharbi
%*/
function [solution,cout] = CS_Finale(cc, ii ,t , p ,a)%parameters : path to text file; nbr_of_iteration ; Nbr of CS; Probability a ; alpha 
[solution,cout]=CS(getMatriceDistances(cc),ii, t , p ,a);
end
function [solution,cout] = CS_Algorithme_PVC(matriceDistance , nbrIterations , TaillePopulation , Pa , alpha)
global matCout;
matCout=matriceDistance;
global ns;
ns=size(matriceDistance,1);
global N;
N=TaillePopulation;
MeilleursSolutions=[];
populationInitialleCoucous=[];  
for i=1:N
    sol=randperm(ns);
    sol=[sol sol(1,1)];
    populationInitialleCoucous=[populationInitialleCoucous;sol];    
end
while( nbrIterations>0)    
    poussin=[];
    for i=1:size(populationInitialleCoucous,1)       
        poussin = getPoussin(populationInitialleCoucous(i,:),alpha);        
        if(getVecteurLongeur(poussin) < getVecteurLongeur(populationInitialleCoucous(i,:)))
            populationInitialleCoucous(i,:)= poussin;            
        end        
    end        
    ss=selectionMeilleurSolution(populationInitialleCoucous);    
    MeilleursSolutions=[MeilleursSolutions;ss];    
    for i=1:size(populationInitialleCoucous,1)
        if(rand(1)<=Pa)             
            tomporaire=populationInitialleCoucous(i,:);
            indice1=randi([2,ns],1,1);
            indice2=randi([2,ns],1,1);
            tmp=tomporaire(1,indice1);
            tomporaire(1,indice1)=tomporaire(1,indice2);
            tomporaire(1,indice2)=tmp;                
            if(getVecteurLongeur(tomporaire) < getVecteurLongeur(populationInitialleCoucous(i,:)))
                populationInitialleCoucous(i,:)= tomporaire;
            end
        end       
    end        
     ss=selectionMeilleurSolution(populationInitialleCoucous);    
     MeilleursSolutions=[MeilleursSolutions;ss];    
     nbrIterations= nbrIterations - 1; 
end
[solution,cout]=solutionFinale(MeilleursSolutions);
end
function Meilleur=selectionMeilleurSolution(population)
vectDistt=getVecteurLongeur(population)';    
 [vv,ind]=min(vectDistt);
Meilleur = population(ind,:);
end
function [Meilleur,cout]=solutionFinale(population)
vectDistt=getVecteurLongeur(population)';    
 [vv,ind]=min(vectDistt);
Meilleur = population(ind,:);
cout=vv;
end
function vecteurX = getPoussin(vect,alpha)
global ns;
levyVect = DistributionLevy(1,ns+1,1);
newVect = vect + alpha * levyVect ; 
vecteurX= getSolutionAcceptable(newVect); 
end
function Xvalide = getSolutionAcceptable(vecteurXreelPo)
global ns;
newVect=abs(vecteurXreelPo); 
        villes=1:ns;
        vv=[];       
        for i=1:size(newVect,2)-1 
            elem=round(newVect(i));
            diff=abs(elem - villes);
            [val,indice]=min(diff);
           vv=[vv villes(1,indice)];            
            villes(indice)=[];            
        end   
        vv=[vv,vv(1,1)];       
Xvalide=vv;
end
function vectDistancesGlobales=getVecteurLongeur(matResultatt) 
global matCout;
vecttt=[];
for elem=matResultatt'
    dist=0;
    for i=1:size(elem,1)-1
        dist = dist+matCout(elem(i,1),elem(i+1,1));
        
    end
    vecttt=[vecttt;dist];
end
vectDistancesGlobales=vecttt;
end
function z = DistributionLevy(n,m,beta)
    num = gamma(1+beta)*sin(pi*beta/2);     
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2); 
    sigma_u = (num/den)^(1/beta);  
    u = random('Normal',0,sigma_u^2,n,m);     
    v = random('Normal',0,1,n,m);
    z = u./(abs(v).^(1/beta));    
end
function matriceDistances =getMatriceDistances(path)
idF=fopen(path);
first=fscanf(idF,'%f',[1 1])
ns=first(1,1);
lecteur=fscanf(idF,'%f',[2,ns]);
coordonne=lecteur'
global coordonner;
coordonner = coordonne;
mat=zeros(ns);
for i=1:ns   
    for j=1:ns        
       mat(i,j)=sqrt(((coordonne(i,1)-coordonne(j,1))^2)+((coordonne(i,2)-coordonne(j,2))^2));       
    end   
end
fclose(idF);
matriceDistances=mat;
end
