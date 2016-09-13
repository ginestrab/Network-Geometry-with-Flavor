function [a,kn,kt,SC] = NGF_d3(N,s,beta,figure_l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code, please cite 
% G. Bianconi and C. Rahmede 
% "Network geometry with flavour: from complexity to quantum geometry"
%Physical Review E 93, 032315 (2016). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code that generates NGF in dimension d=3 and flavour s=-1,0,1.

% a adjacency matrix
% kn vector of  generalized degrees k_{3,0}  of the nodes
% kl vector of  generalized degrees k_{3,1} of links  
% kt vector of  generalized degrees k_{3,2} of triangles  
% SC cell arry of the simplices

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

nt=0;   
SC=cell(1,4);
SC{1}=([1:N])';
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
r(1,:)=[1,1,1]/sqrt(3);
r(2,:)=[-1,-1,1]/sqrt(3);
r(3,:)=[-1,1,-1]/sqrt(3);
r(4,:)=[1,-1,-1]/sqrt(3);
for i=1:4,
    d_N(i)=1;
    D(i)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Initial condition: at time t=1 a single tedrahedron (1,2,3,4)
ntt=1;
SC{4}(ntt,1)=1;
SC{4}(ntt,2)=2;
SC{4}(ntt,3)=3;
SC{4}(ntt,4)=4;
for i1=1:4,
    for i2=(i1+1):4,        
        a(i1,i2)=1;
        a(i2,i1)=1; 
       for i3=(i2+1):4,           
           nt=nt+1;
           SC{3}(nt,1)=i1;
            SC{3}(nt,2)=i2;
             SC{3}(nt,3)=i3;
           tri(nt,1)=i1;
           tri(nt,2)=i2;
           tri(nt,3)=i3;
           at(nt)=exp(-beta*(epsilon(i1)+epsilon(i2)+epsilon(i3)));
           a_occ(nt)=1;
           a_occ3(nt)=1;
       end
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At each time t=in-3 we attach a new tetrahedron

for in=4+1:N,
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
 
    l(1)=tri(nj,1);
    l(2)=tri(nj,2);
    l(3)=tri(nj,3);
    d_N(in)=min([d_N(l(1)),d_N(l(2)),d_N(l(3))])+1;
    D(in)=max(d_N);
    r(in,:)=r(l(1),:)+r(l(2),:)+r(l(3),:);
    r(in,:)=r(in,:)/sqrt(r(in,1)^2+r(in,2)^2+r(in,3)^2);
     a_occ(nj)=a_occ(nj)+s;
     a_occ3(nj)=a_occ3(nj)+1;
    %Add the tethaedron
    for n=1:3,
    a(in,l(n))=1;
    a(l(n),in)=1;  
    end
    for n1=1:3,
        for n2=n1+1:3,
            a(l(n1),l(n2))=a(l(n1),l(n2))+1;
            a(l(n2),l(n1))=a(l(n2),l(n1))+1;
        end
    end
    ntt=ntt+1;
    SC{4}(ntt,4)=in;
    for n=1:3,
          SC{4}(ntt,n)=l(n);
        for n2=n+1:3,           
            nt=nt+1;
            tri(nt,1)=l(n);
            tri(nt,2)=l(n2);
            tri(nt,3)=in;
            SC{3}(nt,1)=l(n);
            SC{3}(nt,2)=l(n2);
             SC{3}(nt,3)=in;
            at(nt)=exp(-beta*(epsilon(l(n))+epsilon(l(n2))+epsilon(in)));
            a_occ(nt)=1;
            a_occ3(nt)=1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized degrees
kn=sum(a>0);

[I2,J2,A2]=find(tril(a));
for n=1:numel(A2),
    SC{2}(n,1)=I2(n);
    SC{2}(n,2)=J2(n);
end
kl=A2;
kt=a_occ3;
a=a>0;
N2=numel(kn);


kn=kn-2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Print network file
[I2,J2,A2]=find(tril(a>0));
if (figure_l==1)
   filename=sprintf('NGF_edgelist_d3_s%d.edges',s);
   fid=fopen(filename,'w');
   for it=1:numel(A2),
       fprintf(fid, ' %d %d  \n', I2(it), J2(it));
   end
   fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
