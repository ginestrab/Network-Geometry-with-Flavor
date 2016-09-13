function [a,kn] = NGF_d1(N,s,beta,figure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code, please cite 
% G. Bianconi and C. Rahmede 
% "Network geometry with flavour: from complexity to quantum geometry"
%Physical Review E 93, 032315 (2016).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code that generates NGF in dimension d=2 and flavour s=-1,0,1.

% a adjacency matrix
% kn vector of  generalized degrees k_{1,0} (the degree) of the nodes

% This code uses 
% N maximal number of nodes in the NGF
% Flavour of the NGF  s=-1,0,1
% Inverse temperature: beta>0 or beta=0
% figure=1 will print the edge list of the network in file 
% "NGF_edgelist_d1_s%d.edges"
% figure=0 will not print the edge list of the network
% energy of the nodes epsilon is uniform from 0-9

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
a=sparse(N,N);
a_occ=zeros(1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign energies to the nodes
% If using Poisson and power-law you must define
% the parameters mu, or kappa
% Examples:
% mu=10;
% kappa=1;
for i=1:N 
epsilon(i)=floor(10*rand(1));
% Alternative energy distributions
% epsilon(i)=random('Poisson',mu); 
% poisson distribution with average mu
% epsilon(i)=rand(1)^(1/(kappa+1)); 
% power-law distribution with exponent kappa

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Initial condition: at time t=1 a single link (1,2)

a(1,2)=exp(-beta*(epsilon(1)+epsilon(2)));
a(2,1)=exp(-beta*(epsilon(1)+epsilon(2)));
k(1)=1;  
k(2)=1;  
a_occ(1)=1;
a_occ(2)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Addition of new links at time t=in-1 the node in is added to the
% network geometry with flavour

for in=2+1:N,
% Choose the node to which attach a new link

V=exp(-beta*epsilon).*a_occ;
norm=sum(V);   
x=rand(1)*norm;
if (norm>0)
for nj1=1:in-1,
x=x-V(nj1);
if x<0,
j=nj1;
break;
end
end
end
% Attach the new link between node in and node j
a(in,j)=exp(-beta*epsilon(in)-beta*epsilon(j));
a(j,in)=a(in,j);
a_occ(in)=1;
a_occ(j)=a_occ(j)+s;	  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized degree (degree) of the nodes
kn=sum(a>0);
a=a>0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print network file
if figure==1
[I,J,A]=find(tril(a));
filename=sprintf('NGF_edgelist_d1_s%d.edges',s);
fid=fopen(filename,'w');
for it=1:max(size(A)),
fprintf(fid, '%d  %d  \n', I(it), J(it));
end
fclose(fid);
end

end
