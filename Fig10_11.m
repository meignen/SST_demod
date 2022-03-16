function Fig11()
% Tests the computation of basin of attraction on the real bat signal
% located at ../data/batsig.txt.


load -ascii batsig.txt
s = batsig;
s = s(150:end)';
s = hilbert(s);

N = length(s);

t =(0:N-1)/N;
gamma = 1e-6;
sigma_opt = 0.1;

clwin = 5;
nr = 3;
lambda = 0;
beta   = 0;

%%computation of transforms with different frequency resolutions
if (N == 2^(floor(log2(N))))
 Ntilde = N;
else
 Ntilde = 2^(floor(log2(N))+1);
end

Nfft1 = Ntilde;
%[STFT,SST,VSST_s] = sst2_new(s,1/sigma_opt^2/N,Nfft1,gamma);
[STFT,SST,VSST_s] = sst2(s,sigma_opt,Nfft1,gamma);

Nfft2 = 2*Ntilde;
%[STFT1,SST1,VSST1_s] = sst2_new(s,1/sigma_opt^2/N,Nfft2,gamma);
[STFT1,SST1,VSST1_s] = sst2(s,sigma_opt,Nfft2,gamma);

Nfft3 = 4*Ntilde;
%[STFT2,SST2,VSST2_s] = sst2_new(s,1/sigma_opt^2/N,Nfft3,gamma);
[STFT2,SST2,VSST2_s] = sst2(s,sigma_opt,Nfft3,gamma);

Nfft4 = 8*Ntilde;
%[STFT3,SST3,VSST3_s] = sst2_new(s,1/sigma_opt^2/N,Nfft4,gamma);
[STFT3,SST3,VSST3_s] = sst2(s,sigma_opt,Nfft4,gamma);

[Cs2, Es] = exridge_mult(VSST_s,nr,lambda,beta,clwin);

[sp1_s,sp2_s,sp3_s,sp4_s,integ1,integ2,integ3,integ4] = ...
 demod_multi(s,VSST_s,VSST1_s,VSST2_s,VSST3_s,Nfft1,Nfft2,Nfft3,Nfft4,t,N,nr,clwin);



d = 0:10;
sign_direct = zeros(length(d),3,N);
sign_demod  = zeros(length(d),3,3,N);

snr1 = [];
snr2 = [];
snr3 = [];
snr4 = [];

s1 = real(s);
X = zeros(1,N);
X1 = zeros(1,N);
X2 = zeros(1,N);
X3 = zeros(1,N);
 
for k = 1:length(d)
 for p = 1:nr
   p
   sign_direct(k,p,:) = 1/Nfft1*real(recmodes(VSST_s,Cs2(p,:),d(k)));   
   %[STFT1,SST1,VSST_sd] = sst2_new(sp1_s(p,:),1/sigma_opt^2/N,Nfft1,gamma);
   [STFT1,SST1,VSST_sd] = sst2(sp1_s(p,:),sigma_opt,Nfft1,gamma);
  
   VSST_sd_int        = zeros(size(VSST_sd));
   A = size(VSST_sd);
   VSST_sd_int(max(1,floor(Nfft1/N*(100-clwin))):min(A(1),floor(Nfft1/N*(100+clwin))),:) = ...
       VSST_sd(max(1,floor(Nfft1/N*(100-clwin))):min(A(1),floor(Nfft1/N*(100+clwin))),:);
   
   [C, Es] = exridge(VSST_sd_int,lambda,beta,clwin);  
   imf                = 1/Nfft1*recmodes(VSST_sd,C,d(k));
   sign_demod(k,p,1,:)= real(imf.*exp(2*1i*pi*(integ1(p,:)-100.*t)));
   
   %[STFT1,SST1,VSST_sd] = sst2_new(sp3_s(p,:),1/sigma_opt^2/N,Nfft1,gamma);
   [STFT1,SST1,VSST_sd] = sst2(sp3_s(p,:),sigma_opt,Nfft1,gamma);
   
   VSST_sd_int        = zeros(size(VSST_sd));
   VSST_sd_int(max(1,floor(Nfft1/N*(100-clwin))):min(A(1),floor(Nfft1/N*(100+clwin))),:) = ...
       VSST_sd(max(1,floor(Nfft1/N*(100-clwin))):min(A(1),floor(Nfft1/N*(100+clwin))),:);
   [C, Es] = exridge(VSST_sd_int,lambda,beta,clwin);  
   imf                    = 1/Nfft1*recmodes(VSST_sd,C,d(k));
   sign_demod(k,p,2,:)    = real(imf.*exp(2*1i*pi*(integ3(p,:)-100.*t)));
   
   %[STFT1,SST1,VSST_sd] = sst2_new(sp4_s(p,:),1/sigma_opt^2/N,Nfft1,gamma);
   [STFT1,SST1,VSST_sd] = sst2(sp4_s(p,:),sigma_opt,Nfft1,gamma);
   
   VSST_sd_int        = zeros(size(VSST_sd));
   VSST_sd_int(max(1,floor(Nfft1/N*(100-clwin))):min(A(1),floor(Nfft1/N*(100+clwin))),:) = ...
       VSST_sd(max(1,floor(Nfft1/N*(100-clwin))):min(A(1),floor(Nfft1/N*(100+clwin))),:);
   [C, Es]                = exridge(VSST_sd_int,lambda,beta,clwin);  
   imf                    = 1/Nfft1*recmodes(VSST_sd,C,d(k));
   sign_demod(k,p,3,:)    = real(imf.*exp(2*1i*pi*(integ4(p,:)-100.*t)));
 end
 X(:)  = sum(sign_direct(k,:,:));
 X1(:) = sum(sign_demod(k,:,1,:));
 X2(:) = sum(sign_demod(k,:,2,:));
 X3(:) = sum(sign_demod(k,:,3,:));
 
 snr1  = [snr1 snr(s1,X-s1)];
 snr2  = [snr2 snr(s1,X1-s1)];
 snr3  = [snr3 snr(s1,X2-s1)];
 snr4  = [snr4 snr(s1,X3-s1)];
end

 figure()
 fs = N/Nfft1*(0:Nfft1/2-1);
 imagesc(t,fs,abs(VSST_s(1:Nfft1/2,:)));
 colormap(1-gray)
 set(gca,'ydir','normal');
 hold on;
 plot(t,N/Nfft1*(Cs2(1,:)-1),t,N/Nfft1*(Cs2(2,:)-1),'--',t,N/Nfft1*(Cs2(3,:)-1),'-.')
 hold off;
 
 figure()
 plot(d,snr1,d,snr2,'--',d,snr3,':',d,snr4,'-o');
 
 [imf] = emd(s1);
 figure()
 %[STFT,SST,VSST_s]    = sst2_new(imf(1,:),1/sigma_opt^2/N,Nfft1,gamma);
 [STFT,SST,VSST_s] = sst2(imf(1,:),sigma_opt,Nfft1,gamma);
 imagesc(t,fs,abs(VSST_s(1:Nfft1/2,:)));
 colormap(1-gray)
 set(gca,'ydir','normal');
 
 figure()
 %[STFT,SST,VSST_s]    = sst2_new(imf(2,:),1/sigma_opt^2/N,Nfft1,gamma);
 [STFT,SST,VSST_s] = sst2(imf(2,:),sigma_opt,Nfft1,gamma);
 imagesc(t,fs,abs(VSST_s(1:Nfft1/2,:)));
 colormap(1-gray)
 set(gca,'ydir','normal');
 
 figure()
 %[STFT,SST,VSST_s]    = sst2_new(imf(3,:),1/sigma_opt^2/N,Nfft1,gamma);
 [STFT,SST,VSST_s] = sst2(imf(3,:),sigma_opt,Nfft1,gamma);
 imagesc(t,fs,abs(VSST_s(1:Nfft1/2,:)));
 colormap(1-gray)
 set(gca,'ydir','normal');
 
end 