function [a,SC] = NGF_d4(N,s,beta,figure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code, please cite 
% G. Bianconi and C. Rahmede 
% "Network geometry with flavour: from complexity to quantum geometry"
%Physical Review E 93, 032315 (2016). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code that generates NGF in dimension d=3 and flavour s=-1,0,1.

% a adjacency matrix
%SC cell array of the simplicies

% This code uses 
% N maximal number of nodes in the NGF
% Flavour of the NGF  s=-1,0,1
% Inverse temperature: beta>0 or beta=0
% figure=1 will print the edge list of the network in file 
%"NGF_edgelist_d3_s%d.edges"
% figure=0 will not print the edge list of the network
% energy of the nodes epsilon is uniform from 0-9

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
a=sparse(N,N);
SC=cell(1,5);
SC{1}=([1:N])';
nt=0;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign energies to the nodes
% If using Poisson and power-law you must define
% the parameters mu, or kappa
% Examples:
% mu=10;
% kappa=1;
for i=1:N,
     epsilon(i)=floor(10*rand(1));
% Alternative energy distributions
    %epsilon(i)=random('Poisson',mu); 
    %poisson distribution with average mu
    %epsilon(i)=rand(1)^(1/(kappa+1)); 
    %power-law distribution with exponent kappa
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Initial condition: at time t=1 a single tedrahedron (1,2,3,4)
L5=1;
SC{5}(L5,1)=1;
SC{5}(L5,2)=2;
SC{5}(L5,3)=3;
SC{5}(L5,4)=4;
SC{5}(L5,5)=5;
L2=0;
L3=0;
for i1=1:5,
    for i2=(i1+1):5,        
        a(i1,i2)=1;
        a(i2,i1)=1;
        L2=L2+1;
        SC{2}(L2,1)=i1;
        SC{2}(L2,2)=i2;
       for i3=(i2+1):5,
             L3=L3+1;
        SC{3}(L3,1)=i1;
        SC{3}(L3,2)=i2;
        SC{3}(L3,3)=i3;
          for i4=(i3+1):5, 
              nt=nt+1;
              qdr(nt,1)=i1;
              qdr(nt,2)=i2;
              qdr(nt,3)=i3;
              qdr(nt,4)=i4;
              SC{4}(nt,1)=i1;
              SC{4}(nt,2)=i2;
              SC{4}(nt,3)=i3;
              SC{4}(nt,4)=i4;
              at(nt)=exp(-beta*(epsilon(i1)+epsilon(i2)+epsilon(i3)+epsilon(i4)));
              a_occ(nt)=1;
              a_occ3(nt)=1;
          end
       end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At each time t=in-3 we attach a new tetrahedron

for in=5+1:N,
% Choose triangular face to which to attach the new tetrahedron 

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
 
    l(1)=qdr(nj,1);
    l(2)=qdr(nj,2);
    l(3)=qdr(nj,3);
    l(4)=qdr(nj,4);
     a_occ(nj)=a_occ(nj)+s;
     a_occ3(nj)=a_occ3(nj)+1;
    %Add the tethaedron
    L5=L5+1;
     SC{5}(L5,1)=l(1);
        SC{5}(L5,2)=l(2);
      SC{5}(L5,3)=l(3);
        SC{5}(L5,4)=l(4);
          SC{5}(L5,5)=in;

    
    for n=1:4,
        L2=L2+1;
        SC{2}(L2,1)=l(n);
        SC{2}(L2,2)=in;
    a(in,l(n))=1;
    a(l(n),in)=1;  
    end
    for n1=1:4,
        for n2=n1+1:4,
            a(l(n1),l(n2))=a(l(n1),l(n2))+1;
            a(l(n2),l(n1))=a(l(n2),l(n1))+1;
        end
    end
    for n=1:4,
        for n2=n+1:4,
            L3=L3+1;
            SC{3}(L3,1)=l(n);
                SC{3}(L3,2)=l(n2);
                SC{3}(L3,3)=in;
            for n3=n2+1:4,
                nt=nt+1;
                qdr(nt,1)=l(n);
                qdr(nt,2)=l(n2);
                qdr(nt,3)=l(n3);
                qdr(nt,4)=in;
                SC{4}(nt,1)=l(n);
                SC{4}(nt,2)=l(n2);
                SC{4}(nt,3)=l(n3);
                SC{4}(nt,4)=in;
                at(nt)=exp(-beta*(epsilon(l(n))+epsilon(l(n2))+epsilon(l(n3))+epsilon(in)));
                a_occ(nt)=1;
                a_occ3(nt)=1;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized degrees
kn=sum(a>0);
kn=kn-2;
[I2,J2,A2]=find(tril(a));
kl=A2;
kt=a_occ3;
a=a>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Print network file
[I2,J2,A2]=find(tril(a>0));
if (figure==1)
   filename=sprintf('NGF_edgelist_d4_s%d.edges',s);
   fid=fopen(filename,'w');
   for it=1:numel(A2),
       fprintf(fid, ' %d %d  \n', I2(it), J2(it));
   end
   fclose(fid);
end

