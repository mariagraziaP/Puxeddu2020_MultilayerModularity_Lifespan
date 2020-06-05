function [Fmat, nodF, layF] = myFlexibility(part)

% input: part = nodes*layers -> matrix of communities

N = size(part,1);
nL = size(part,2);

Fmat = zeros(N,nL-1);
for idx=1:(nL-1)
    vec = part(:,idx)-part(:,idx+1);
    Fmat(vec~=0,idx) = 1;
    clear vec
end
clear idx

nodF = sum(Fmat,2)/(nL-1);

layF = sum(Fmat,1);