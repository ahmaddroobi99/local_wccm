function [e,kp,count] = spectra_omnidirectional_sum(qhcv)
% display(qhcv);
N = size(qhcv,1);
% display(N);
k = [0:N/2 -N/2+1:-1]';
% disp(k);
[KX, KY] = meshgrid(k,k);
% disp(KX);
% disp(KY);
Knorm = hypot(KX,KY);
% display(Knorm);


e = zeros(N/2+1,1);
% display(e);

count = zeros(N/2+1,1);
% display(count);

for jj = 1:N
    for ii = 1:N
        k = Knorm(ii,jj);
        if ceil(k) <= N/2
            r = k-floor(k);
            % plus 1 because of matlab indexing
            count(floor(k)+1) =  count(floor(k)+1) + 1;
            count(ceil(k)+1) =  count(ceil(k)+1) + 1;
            
            e(floor(k)+1) = e(floor(k)+1) + (1-r)*abs(qhcv(ii,jj));
            e(ceil(k)+1) = e(ceil(k)+1) + r*abs(qhcv(ii,jj));
        end
    end
end
e = 0.5*e;
% display(e);

kp = Knorm(1,1:N/2+1);
% display(kp);


count = 0.5*count;
% display(count);

