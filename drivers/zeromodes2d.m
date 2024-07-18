% function inds = zeromodes2d(N)
% % return a vector of the indices corresponding to the zero modes of a N x N
% % matrix 
% 
% inds = zeros(4,1); 
% inds(1) = sub2ind([N,N],1,1);
% inds(2) = sub2ind([N,N],N/2+1,1);
% inds(3) = sub2ind([N,N],1, N/2+1);
% inds(4) = sub2ind([N,N],N/2+1, N/2+1);
% 
% 


% Function to get zero modes
function inds = zeromodes2d(N)
% Return a vector of the indices corresponding to the zero modes of an N x N matrix 

inds = zeros(4,1); 
inds(1) = sub2ind([N,N],1,1);
inds(2) = sub2ind([N,N],N/2+1,1);
inds(3) = sub2ind([N,N],1, N/2+1);
inds(4) = sub2ind([N,N],N/2+1, N/2+1);
end


