function [xxq, yyq, wq] = vorticity(qh, ntimes)
% compute vorticity from fourier coeffs of voritcity qh
% ntimes: number of times to interpolate data

if nargin < 2
    ntimes = 4;
end

L = 2*pi;
N = size(qh,1);

xq = 0:1/(2^ntimes):N-1; xq = L/N * xq;
yq = 0:1/(2^ntimes):N-1; yq = L/N * yq;
[xxq, yyq] = meshgrid(xq,yq);

kx = [0:N/2-1 0 -N/2+1:-1]; % x wave vector
ky = [0:N/2-1 0 -N/2+1:-1]; % y wave vector
[kxx, kyy] = meshgrid(kx, ky);
 
qh(zeromodes2d(N)) = 0;
w = ifft2(qh*N*N);

wq = interp2(w, ntimes, 'cubic');
wq = real(wq); % small numerical imaginary parts may accumulate
end
