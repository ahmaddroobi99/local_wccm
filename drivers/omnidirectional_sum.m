function [e,kp,count] = omnidirectional_sum(qhcv)

if size(qhcv,1) == 1 || size(qhcv,2) == 1
    qhcv = reshape(qhcv,[1,1]*sqrt(size(qhcv,1)));
end

N = size(qhcv,1);
k = [0:N/2 -N/2+1:-1]';
[KX, KY] = meshgrid(k,k);
Knorm = hypot(KX,KY);

e = zeros(N/2+1,1);
count = zeros(N/2+1,1);

for jj = 1:N
    for ii = 1:N
        k = Knorm(ii,jj);
        if ceil(k) <= N/2
            r = k-floor(k);
            % plus 1 because of matlab indexing
            count(floor(k)+1) =  count(floor(k)+1) + (1-r);
            count(ceil(k)+1) =  count(ceil(k)+1) + r;
            
            e(floor(k)+1) = e(floor(k)+1) + (1-r)*abs(qhcv(ii,jj));
            e(ceil(k)+1) = e(ceil(k)+1) + r*abs(qhcv(ii,jj));
        end
    end
end

e = e./count;
kp = Knorm(1,1:N/2+1);