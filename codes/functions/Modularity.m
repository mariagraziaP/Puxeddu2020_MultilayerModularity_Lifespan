function Q = Modularity(adj,com,gamma,type)

% adj = adjacency matrix
% com = community membership
% type = 'dir' (directed network) 'und' (undirected network)

N = length(adj);

switch type
    case 'dir'
        Ki=sum(adj,1);                            %in-degree
        Ko=sum(adj,2);                            %out-degree
        m=sum(Ki);                           	%number of edges
        b=adj-gamma*(Ko*Ki).'/m;
        B=b+b.';                            	%directed modularity matrix
        s=com(:,ones(1,N));                      %compute modularity
        Q=~(s-s.').*B/(2*m);
        Q=sum(Q(:));
    case 'und'
        K=sum(adj);                               %degree
        m=sum(K);                               %number of edges (each undirected edge is counted twice)
        B=adj-gamma*(K.'*K)/m;                    %modularity matrix
        s=com(:,ones(1,N));                      %compute modularity
        Q=~(s-s.').*B/m;
        Q=sum(Q(:));
end





