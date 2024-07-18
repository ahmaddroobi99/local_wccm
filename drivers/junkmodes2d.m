function inds = junkmodes2d(N,enforcezero)
%-------------------------------------------------------------------------%
% return a vector of the indices corresponding to the junk modes of a N x N
% Fourier matrix in the full plane, i.e. the row and column with the
% highest N/2+1 modes that are zero
%
% Warning: enforce 0 for the zero mode!
%-------------------------------------------------------------------------%
if nargin <2
    enforcezero = true;
end

picker = reshape(1:N*N, N, N);
inds = [picker(:,N/2+1); picker(N/2+1,1:N/2)'; picker(N/2+1,N/2+2:N)']; % make sure we dont double include the corner

if enforcezero
    inds = [1; inds];
end
