function [contrQ,listComm] = contribution2modularity(adj,com,gamma)

% NB: only for undirected networks

N = size(adj,1);
K=sum(adj);                                 %degree
m=sum(K);                                   %number of edges (each undirected edge is counted twice)
B=adj-gamma*(K.'*K)/m;                    	%modularity matrix
s=com(:,ones(1,N));                         %compute modularity
Q=~(s-s.').*B/m;

lab = unique(com);
nc = length(lab);
for n=1:nc
    nodes = find(com==lab(n));
    subQ = Q(nodes,nodes);
    contrQ(n) = sum(subQ(:));
    listComm{n,1} = lab(n);
    listComm{n,2} = nodes;
    
    clear nodes subQ
end
clear nc

