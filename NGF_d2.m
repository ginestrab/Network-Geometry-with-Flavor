

function [a,kn,kl,SC] = NGF_d2(N,s,beta,figure_l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code, please cite 
% G. Bianconi and C. Rahmede 
% "Network geometry with flavour: from complexity to quantum geometry"
%Physical Review E 93, 032315 (2016). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code that generates NGF in dimension d=2 and flavour s=-1,0,1.

% a adjacency matrix
% kn vector of  generalized degrees k_{2,0}  of the nodes
% kl vector of  generalized degrees k_{2,1} of links  
% SC cell array of the simplicies

% This code uses 
% N maximal number of nodes in the NGF
% Flavour of the NGF  s=-1,0,1
% Inverse temperature: beta>0 or beta=0
% figure=1 will print the edge list of the network in file 
% "NGF_edgelist_d2_s%d.edges"
% figure=0 will not print the edge list of the network
% energy of the nodes epsilon is uniform from 0-9

SC=cell(1,3);
SC{1}=([1:N])';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inizialization

a=sparse(N,N);
a_occ=sparse(N,N);
a_occ2=sparse(N,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign energies to the nodes
% If using Poisson and power-law you must define 
% the parameters mu, or kappa
% Examples:
% mu=10;
% kappa=1
for i=1:N
     epsilon(i)=floor(10*rand(1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition at time t=1 including a single triangle between nodes
% 1,2,3
L=0;
for i1=1:3,
    
    for i2=(i1+1):3,     
        L=L+1;
        a(i1,i2)=exp(-beta*(epsilon(i1)+epsilon(i2)));
        a(i2,i1)=exp(-beta*(epsilon(i1)+epsilon(i2)));
        a_occ(i1,i2)=1;  
        a_occ(i2,i1)=1;  
        a_occ2(i1,i2)=1;
        a_occ2(i2,i1)=1;
    end
end
nt=1;

    SC{3}(nt,1)=1;
     SC{3}(nt,2)=2;
      SC{3}(nt,3)=3;

for i=1:3,
    r(i,:)=[cos(2*pi*(i-1)/3),sin(2*pi*(i-1)/3)];
    d_N(i)=1;
    D(i)=1;
    theta(i)=1/3*(i-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At each time t=in-2 we attach a new triangle

for in=(3+1):N,
    % Choose edge (l1,l2) to which we will attach the new triangle
    
    [I,J,V]=find(tril(a.*(a_occ)));
   
    norm=sum(V);   
    x=rand(1)*norm;
    if (norm>0)
	    for nj1=1:numel(V),
		    x=x-V(nj1);
            if x<0,
                nj=nj1;
                break;
            end
        end
	    l1=I(nj);
	    l2=J(nj);
        
        d_N(in)=min(d_N(l1),d_N(l2))+1;
        D(in)=max(d_N);
        r(in,:)=r(l1,:)+r(l2,:);
        r(in,:)=r(in,:)/sqrt(r(in,1)^2+r(in,2)^2);
       theta(in)=min(theta(l1),theta(l2))+1/2*abs(theta(l1)-theta(l2));
        if((theta(l1)==0)&&(theta(l2)>=2/3)),
            theta(in)=theta(l2)+1/2*(1-theta(l2));
        end
        if((theta(l2)==0)&&(theta(l1)>=2/3)),
            theta(in)=theta(l1)+1/2*(1-theta(l1));
        end
        
   
	     a_occ(l1,l2)=a_occ(l1,l2)+s;
	     a_occ(l2,l1)=a_occ(l2,l1)+s;
	     a_occ2(l1,l2)=a_occ2(l1,l2)+1;
         a_occ2(l2,l1)=a_occ2(l2,l1)+1;
      
         nt=nt+1;
    SC{3}(nt,1)=in;
     SC{3}(nt,2)=l1;
      SC{3}(nt,3)=l2;
      
         % Attach the new node in to the node l1;
         L=L+1;
	     a(in,l1)=exp(-beta*(epsilon(l1)+epsilon(in)));
         a(l1,in)=exp(-beta*(epsilon(l1)+epsilon(in)));
         a_occ(in,l1)=1;
         a_occ(l1,in)=1;
         a_occ2(in,l1)=1;
         a_occ2(l1,in)=1;
        
         % Attach the new node in to the node l2;
         L=L+1;
         a(in,l2)=exp(-beta*(epsilon(l2)+epsilon(in)));
         a(l2,in)=exp(-beta*(epsilon(l2)+epsilon(in)));
         a_occ(in,l2)=1;
         a_occ(l2,in)=1;
         a_occ2(in,l2)=1;
         a_occ2(l2,in)=1;
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized degrees 
k=sum(a>0);

[I,J,kl]=find(tril(a_occ2));
kn=k-ones(size(k));
a=a>0;

[I,J,V]=find(tril(a>0));
for n=1:numel(V),
SC{2}(n,1)=I(n);
SC{2}(n,2)=J(n);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Print network file

if figure_l==1
    [I,J,A]=find(tril(a));
    filename=sprintf('NGF_edgelist_d2_s%d.edges',s);
    fid=fopen(filename,'w');
    for it=1:max(size(A)),
        fprintf(fid, '%d  %d  \n', I(it), J(it));
    end
    fclose(fid);
end


end



