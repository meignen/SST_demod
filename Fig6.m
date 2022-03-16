function [En] = Fig7(cas)

 %% Signal
 %N = 2048;
 N = 1024;
 t = (0:N-1)/N;

 %test signal 1
 if (cas == 2)
 %test signal 1
   a = 2;
   s1 = a.*exp(2*pi*1i*(330*t+16*cos(3*pi*t)));
   s2 = a.*exp(2*pi*1i*(190*t+9*cos(3*pi*t)));
   s3 = a.*exp(2*pi*1i*(40*t));
   s  = s1+s2+s3;
   clwin = 10;
   nc = 3;
 else
 %test signal 2
  a  = 2;
  s1 = a.*exp(2*pi*1i*(250*t+50*t.^3));
  s2 = a.*exp(2*pi*1i*(130*t+100*t.^2));
  s3 = a.*exp(2*pi*1i*(90*t+0.2*cos(3*pi*t)));%low frequency component
  s  = s1+s2+s3;
  clwin = 30;
  nc = 3;
 end

 %% Parameters
 gamma = 0.01;
 sigma = 0.005:0.01:0.1;

 En  = zeros(4,length(sigma));

 lambda = 0;
 beta   = 0;
 
 for p = 1:length(sigma)
  Nfft1 = N;
  [STFT,SST,VSST_s] = sst2(s,sigma(p),Nfft1,gamma);
  [Cs2, Es] = exridge_mult(VSST_s,nc,lambda,beta,clwin);
  En(1,p) = sum(Es)/sum((abs(VSST_s(:))).^2);
 end

 %noisy case
 SNR = 5:-5:-5;
 for k = 1:length(SNR)
  k
  En_int  = zeros(5,length(sigma));
  for q = 1:10
   q
   s0 = add_awgn_noise(s,SNR(k));
   for p = 1:length(sigma)
    Nfft1 = N;
    [STFT,SST,VSST_s] = sst2(s0,sigma(p),Nfft1,gamma);
    [Cs2, Es] = exridge_mult(VSST_s,nc,lambda,beta,clwin);
    En_int(q,p) = sum(Es)/sum((abs(VSST_s(:))).^2);
   end
  end 
  En(k+1,:)     = mean(En_int);
 end
 figure()
 plot(sigma,En(1,:),sigma,En(2,:),'--',sigma,En(3,:),':',sigma,En(4,:),'-o');
end