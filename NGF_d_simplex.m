function [a,k] = NGF_d_simplex(N,s,beta,d,figure_l)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  If you use this code, please cite:
%  G. Bianconi and C. Rahmede
%  "Network geometry with flavour: from complexity to quantum geometry"
%  Physical Review E 93, 032315 (2016).
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code that generates NGF in dimension d and flavour s=-1,0,1 costructed
% with d-dimensional SIMPLICES

% a adjacency matrix
% k vector of  degrees   of the nodes


% This code uses
% N maximal number of nodes in the NGF
% Flavour of the NGF  s=-1,0,1
% Inverse temperature: beta>0 or beta=0
% Dimension d with d>1
% figure_l=1 will print the edge list of the network in file
%"NGF_edgelist_d%d_s%d.edges"
% figure_l=0 will not print the edge list of the network
% energy of the nodes epsilon is uniform from 0-9% Code that generates NGF in dimension d and flavour s=-1,0,1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign energies to the nodes
% If using Poisson and power-law you must define
% the parameters mu, or kappa
% Examples:
% mu=10;
% kappa=1;

     epsilon=floor(10*rand(N,1));
% Alternative energy distributions
    %epsilon(i)=random('Poisson',mu); 
    %poisson distribution with average mu
    %epsilon(i)=rand(1)^(1/(kappa+1)); 
    %power-law distribution with exponent kappa
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Initial condition: at time t=1 a single d-dimensional hypercube (1,2,3,4)


    for i = 1:(d+1)
        for j = (i+1):(d+1)
            a(i,j)=1;
            a(j,i)=1;
        end
    end
        
    for nt = 1:(d+1)
        at(nt)=1;
        a_occ(nt)=1;
        a_occ3(nt)=1;
        j=0;
        for i=1:(d+1),
            if(abs(i-nt)>0)
                j=j+1;
                node(nt,j)=i;
                at(nt)=at(nt)*exp(-beta*epsilon(i));
            end
        end
    end


      it=d+1;  
      
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


while it<N,
    it=it+1;

    [I,J,V]=find(at.*a_occ);

    norm=sum(V);
    x=rand(1)*norm;
    for nj1=1:numel(V),
            x=x-V(nj1);
         if x<0,
             nj=J(nj1);
             break;
         end
    end    


    
     a_occ(nj)=a_occ(nj)+s;
     a_occ3(nj)=a_occ3(nj)+1;
     
     for n1=1:d,
         j=node(nj,n1);
         a(it,j)=1;
         a(j,it)=1;
     end
    for n1=1:d,
        nt=nt+1;
        at(nt)=1;
        a_occ(nt)=1;
        a_occ3(nt)=1;
        node(nt,1)=it;
        at(nt)=1;
        j=1;
        for n2=1:d,
            
            if(abs(n2-n1)>0)
                j=j+1;
             node(nt,j)=node(nj,n2);
           
             at(nt)=at(nt)*exp(-beta*epsilon(node(nj,n2)));
             a(it,node(nj,n2))=1;
             a(node(nj,n2),it)=1;
            end
        end
    end
     

end
    
    k=sum(a>0)-(d-1);

 

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% % Print network file
 [I2,J2,A2]=find(tril(a>0));
 if (figure_l==1)
    filename=sprintf('NGF_edgelist_d%d_s%d.edges',d,s);
    fid=fopen(filename,'w');
    for it=1:numel(A2),
        fprintf(fid, ' %d %d  \n', I2(it), J2(it));
    end
    fclose(fid);
 end
end
