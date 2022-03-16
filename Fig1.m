function Fig1(cas,jump) 
%% Signal
N = 1024;
t = (0:N-1)/N;

a    = 2;
if cas == 1
 %linear chirp
 s = a.*exp(2*pi*1i*(130*t+100*t.^2));
 phi1 = 130+200*t;
elseif cas == 2
 %polynomial chirp    
 s   = a.*exp(2*pi*1i*(250*t+50*t.^3));
 phi1 = 250+150*t.^2;
else
 %sinusoidal wave   
 s = a.*exp(2*pi*1i*(330*t+16*cos(3*pi*t)));
 phi1 = 330-16*3*pi*sin(3*pi*t);
end

%% Parameters
 gamma = 0.01;
 sigma = 0.04;

 % Choice of time and frequency bins
 
 %% Compute TF transforms
 Nfft1 = N;
[STFT,SST,VSST]      = sst2(s,sigma,Nfft1,gamma);
 Nfft2 = 2*N;
[STFT1,SST1,VSST1]   = sst2(s,sigma,Nfft2,gamma);
 Nfft3 = 4*N;
 [STFT2,SST2,VSST2]   = sst2(s,sigma,Nfft3,gamma);
 Nfft4 = 8*N;
 [STFT3,SST3,VSST3]   = sst2(s,sigma,Nfft4,gamma);
 Nfft5 = 16*N;
 [STFT4,SST4,VSST4]   = sst2(s,sigma,Nfft5,gamma);
 
 t  = (0:N-1)/N;
 fs = N/Nfft1*(0:Nfft1/2-1);
 imagesc(t,fs,abs(STFT(1:Nfft1/2,:)));
 set(gca,'YDir','normal')
 
 %Computation of the ridges with no regularization, effect of zero padding 
 beta   = 0;
 lambda = 0;
 [Cs0, Es] = exridge(STFT,lambda,beta,jump);
 [Cs1,Es]  = exridge(SST,lambda,beta,jump);
 [Cs2, Es] = exridge(VSST,lambda,beta,jump);
 [Cs01,Es] = exridge(STFT1,lambda,beta,jump);    
 [Cs11,Es] = exridge(SST1,lambda,beta,jump);
 [Cs21,Es] = exridge(VSST1,lambda,beta,jump);
 [Cs02,Es] = exridge(STFT2,lambda,beta,jump);    
 [Cs12,Es] = exridge(SST2,lambda,beta,jump);
 [Cs22,Es] = exridge(VSST2,lambda,beta,jump);
 [Cs03,Es] = exridge(STFT3,lambda,beta,jump);    
 [Cs13,Es] = exridge(SST3,lambda,beta,jump);
 [Cs23,Es] = exridge(VSST3,lambda,beta,jump);
 [Cs04,Es] = exridge(STFT4,lambda,beta,jump);    
 [Cs14,Es] = exridge(SST4,lambda,beta,jump);
 [Cs24,Es] = exridge(VSST4,lambda,beta,jump);

 %Computation of mean square error
 MSE = zeros(9,5);
 
 %without noise 
 index = 100:N-100;
 %STFT
 MSE(1,:) = [ sqrt(sum((N/Nfft1*(Cs0(index)-1)-phi1(index)).^2)/(length(index)-1)) sqrt(sum((N/Nfft2*(Cs01(index)-1)-phi1(index)).^2)/(length(index)-1))...
              sqrt(sum((N/Nfft3*(Cs02(index)-1)-phi1(index)).^2)/(length(index)-1)) sqrt(sum((N/Nfft4*(Cs03(index)-1)-phi1(index)).^2)/(length(index)-1))...
              sqrt(sum((N/Nfft5*(Cs04(index)-1)-phi1(index)).^2)/(length(index)-1))]; 
 %FSST
 MSE(2,:) = [ sqrt(sum((N/Nfft1*(Cs1(index)-1)-phi1(index)).^2)/(length(index)-1)) sqrt(sum((N/Nfft2*(Cs11(index)-1)-phi1(index)).^2)/(length(index)-1))...
              sqrt(sum((N/Nfft3*(Cs12(index)-1)-phi1(index)).^2)/(length(index)-1)) sqrt(sum((N/Nfft4*(Cs13(index)-1)-phi1(index)).^2)/(length(index)-1))... 
              sqrt(sum((N/Nfft5*(Cs14(index)-1)-phi1(index)).^2)/(length(index)-1))];
 %VSST 
 MSE(3,:) = [ sqrt(sum((N/Nfft1*(Cs2(index)-1)-phi1(index)).^2)/(length(index)-1)) sqrt(sum((N/Nfft2*(Cs21(index)-1)-phi1(index)).^2)/(length(index)-1))...
              sqrt(sum((N/Nfft3*(Cs22(index)-1)-phi1(index)).^2)/(length(index)-1)) sqrt(sum((N/Nfft4*(Cs23(index)-1)-phi1(index)).^2)/(length(index)-1))... 
              sqrt(sum((N/Nfft5*(Cs24(index)-1)-phi1(index)).^2)/(length(index)-1))]; 

%with noise
p = 4;
for k = 5:-5:-5
 k
 MSE_int = zeros(5,2,5); %we consider 5 different realizations 
 for q =1:5
  q
  s0 = add_awgn_noise(s,k);
  [STFT,SST,VSST]    = sst2(s0,sigma,Nfft1,gamma);
  [STFT1,SST1,VSST1] = sst2(s0,sigma,Nfft2,gamma);
  [STFT2,SST2,VSST2] = sst2(s0,sigma,Nfft3,gamma);
  [STFT3,SST3,VSST3] = sst2(s0,sigma,Nfft4,gamma);
  [STFT4,SST4,VSST4] = sst2(s0,sigma,Nfft5,gamma);
  [Cs0, ~] = exridge(STFT,lambda,beta,jump);     
  [Cs2, ~] = exridge(VSST,lambda,beta,jump);
  [Cs01,~] = exridge(STFT1,lambda,beta,jump);   
  [Cs21,~] = exridge(VSST1,lambda,beta,jump);
  [Cs02,~] = exridge(STFT2,lambda,beta,jump);    
  [Cs22,~] = exridge(VSST2,lambda,beta,jump);
  [Cs03,~] = exridge(STFT3,lambda,beta,jump);    
  [Cs23,~] = exridge(VSST3,lambda,beta,jump);
  [Cs04,~] = exridge(STFT4,lambda,beta,jump);    
  [Cs24,~] = exridge(VSST4,lambda,beta,jump);
  
  %STFT
  MSE_int(q,1,:) = [ sqrt(sum((N/Nfft1*(Cs0(index)-1)-phi1(index)).^2)/(length(index)-1)) sqrt(sum((N/Nfft2*(Cs01(index)-1)-phi1(index)).^2)/(length(index)-1))...
                 sqrt(sum((N/Nfft3*(Cs02(index)-1)-phi1(index)).^2)/(length(index)-1)) sqrt(sum((N/Nfft4*(Cs03(index)-1)-phi1(index)).^2)/(length(index)-1))...
                 sqrt(sum((N/Nfft5*(Cs04(index)-1)-phi1(index)).^2)/(length(index)-1))]; 
  %VSST 
  MSE_int(q,2,:) = [ sqrt(sum((N/Nfft1*(Cs2(index)-1)-phi1(index)).^2)/(length(index)-1)) sqrt(sum((N/Nfft2*(Cs21(index)-1)-phi1(index)).^2)/(length(index)-1))...
                     sqrt(sum((N/Nfft3*(Cs22(index)-1)-phi1(index)).^2)/(length(index)-1)) sqrt(sum((N/Nfft4*(Cs23(index)-1)-phi1(index)).^2)/(length(index)-1))... 
                     sqrt(sum((N/Nfft5*(Cs24(index)-1)-phi1(index)).^2)/(length(index)-1))];
 end
 MSE(p,:) = mean(MSE_int(:,1,:));
 p = p+1;
 MSE(p,:) = mean(MSE_int(:,2,:));
 p = p+1;
end

figure()
vec = [1 2 4 8 16];
plot(vec,MSE(1,:),'--',vec,MSE(2,:),vec,MSE(3,:),':',vec,MSE(4,:),'-*',vec,MSE(5,:),'-o',...
     vec,MSE(6,:),'-d',vec,MSE(7,:),'-s',vec,MSE(8,:),'-<',vec,MSE(9,:),'-x');
end
