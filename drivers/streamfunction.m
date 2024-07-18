function [xxq, yyq, Psiq] = streamfunction(U, ntimes)
% compute streamfunction from U
% ntimes: number of times to interpolate data

if nargin < 2
    ntimes = 4;
end

L = 2*pi;
N = size(U,1);

xq = 0:1/(2^ntimes):N-1; xq = L/N * xq;
yq = 0:1/(2^ntimes):N-1; yq = L/N * yq;
[xxq, yyq] = meshgrid(xq,yq);

kx = [0:N/2-1 0 -N/2+1:-1]; % x wave vector
ky = [0:N/2-1 0 -N/2+1:-1]; % y wave vector
[kxx, kyy] = meshgrid(kx, ky);
knorm = hypot(kxx,kyy);

Psih = U./knorm;
Psih(zeromodes2d(N)) = 0;
Psi = ifft2(Psih*N*N);

Psiq = interp2(Psi, ntimes, 'cubic');
Psiq = real(Psiq); % small numerical imaginary parts may accumulate
end
