function [id, idx,idy] = subsample2d(Nt,Ns)



idx = [1:Ns/2+1 Nt-Ns/2+2:Nt];
idy = [1:Ns/2+1 Nt-Ns/2+2:Nt];

ii = 0;
id = zeros(Ns*Ns,1);
for i = idx
    for j = idy
        ii = ii + 1;
        id(ii) = sub2ind([Nt,Nt],j,i);
    end
end
