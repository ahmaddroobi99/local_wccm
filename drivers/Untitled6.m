clear;

N = 4;
N2 = N*N;
Nt = 100;

fname = 'test1.bin';



fid = fopen(fname,'w');
for nn = 1:Nt
    x = ones(N2,1)*nn;
    if ~(mod(nn,1))
        fwrite(fid,x,'double');
    end
end
fclose(fid);


fileID = fopen(fname,'r');
for nn = 1:Nt
A = fread(fileID,[N2],'double') ;
A
end

fclose(fileID);

%%
clc
clear;

N = 2;
N2 = N*N;
Nt = 2;

fname = 'test1.bin';

fid = fopen(fname,'w');
for nn = 1:Nt
    x = ones(N2,1)*nn + 1i*1.1*ones(N2,1)*nn;
    disp(x)
    if ~(mod(nn,1))
        fwrite(fid,[real(x) imag(x)],'double');
    end
end
fclose(fid);
fprintf('\n\n')

fileID = fopen(fname,'r');
for nn = 1:Nt
Ax = fread(fileID, [2*N2], 'double') ;

disp(Ax(1:N2)+ 1i*Ax(N2+1:2*N2));
fprintf('\n')

end

fclose(fileID);








%%
% clear all;

datafolder = fullfile('data','baroRK4_N256tau0_beta0kf50ef20000000000sf1af0_d59.8td10_sinset-02-26-18-113312');

fnameqh = fullfile(datafolder, 'qh.bin');
fnamet = fullfile(datafolder, 't.bin');

dataext = load(fullfile(datafolder, 'diagnostics.mat'));
N =  dataext.params.N;
Ntotal = N*N;

fidqh = fopen(fnameqh,'r');
fidt = fopen(fnamet,'r');

t = fread(fidt,'double');

k = [0:N/2 -N/2+1:-1]';
[KX, KY] = meshgrid(k,k);
dX = 1i*repmat(k', [N 1]);
dY = 1i*repmat(k, [1 N]);
Laplacian = dX.^2 + dY.^2;
InvBT = 1./Laplacian; InvBT(1,1) = 0;

for kk = 1:1:length(t)
    qq = fread(fidqh,2*Ntotal,'double');
    
    qq = complex(qq(1:Ntotal), qq(Ntotal+1:2*Ntotal));
    
    if ~mod(kk,5)   
        qh = reshape(qq,N,N);
        psi_hat = InvBT.*qh;
        pcolor(real(ifft2(psi_hat)));
        shading flat;
        drawnow limitrate
    end
    
end

% ee = zeros(N/2+1,length(t));
% for kk = 1:length(t)
%     qq = fread(fidqh,2*Ntotal,'double');
%     qq = complex(qq(1:Ntotal), qq(Ntotal+1:2*Ntotal));
%     qh = reshape(qq,N,N);
%     psi_hat = InvBT.*qh;
% 
%     
%     ENE = zeros(N/2+1,1);
%     for jj = 1:N
%         for ii = 1:N
%             k = sqrt(KX(ii,jj).^2 + KY(ii,jj).^2);
%             if ceil(k) <= N/2
%                 r = k-floor(k);
%                 % plus 1 because of matlab indexing
%                 ENE(floor(k)+1) = ENE(floor(k)+1) + (1-r)*(k^2)*abs(psi_hat(ii,jj))^2;
%                 ENE(ceil(k)+1) = ENE(ceil(k)+1) + r*(k^2)*abs(psi_hat(ii,jj))^2;
%             end
%         end
%     end
%     ee(:,kk) = .5*ENE/(Ntotal^2);

% end

fclose(fidqh);
fclose(fidt);












