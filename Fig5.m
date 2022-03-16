function [Ren] = Fig6(cas1,SNR)
 
 %cas1: noise or not, noise level given by SNR
 
 %% Signal
 N = 1024;
 t = (0:N-1)/N;

 a  = 2;
 s1 = a.*exp(2*pi*1i*(250*t+50*t.^3));
 s2 = a.*exp(2*pi*1i*(130*t+100*t.^2));
 s3 = a.*exp(2*pi*1i*(90*t+0.2*cos(3*pi*t)));%low frequency component
 s  = s1+s2+s3;

 %% Parameters
 gamma = 0.01;
 sigma = 0.005:0.01:0.1;

 if cas1 == 0 
   Ren = zeros(4,length(sigma));
   for p = 1:length(sigma)
    Nfft1 = N;
    [STFT,SST,VSST_s]    = sst2(s,sigma(p),Nfft1,gamma);
 
    Nfft2 = 2*N;
    [STFT1,SST1,VSST1_s] = sst2(s,sigma(p),Nfft2,gamma);   

    Nfft3 = 4*N;
    [STFT2,SST2,VSST2_s] = sst2(s,sigma(p),Nfft3,gamma); 
 
    Nfft4 = 8*N;
    [STFT3,SST3,VSST3_s] = sst2(s,sigma(p),Nfft4,gamma);
    fs = 0:Nfft1/2-1;
    Ren(1,p) = renyi_entropy(abs(VSST_s(1:Nfft1/2,:)),t,fs',3);
    fs = 0:Nfft2/2-1;
    Ren(2,p) = renyi_entropy(abs(VSST1_s(1:Nfft2/2,:)),t,fs',3);
    fs = 0:Nfft3/2-1; 
    Ren(3,p) = renyi_entropy(abs(VSST2_s(1:Nfft3/2,:)),t,fs',3);
    fs = 0:Nfft4/2-1;
    Ren(4,p) = renyi_entropy(abs(VSST3_s(1:Nfft4/2,:)),t,fs',3);
   end 
 else 
  Ren_bis  = zeros(5,4,length(sigma));
  Ren = zeros(4,length(sigma));
  for q = 1:5
   q   
   s0 = add_awgn_noise(s,SNR);
   for p = 1:length(sigma)
    Nfft1 = N;
   [STFT,SST,VSST_s] = sst2(s0,sigma(p),Nfft1,gamma);

    Nfft2 = 2*N;
    [STFT1,SST1,VSST1_s] = sst2(s0,sigma(p),Nfft2,gamma);   
   
    Nfft3 = 4*N;
    [STFT2,SST2,VSST2_s] = sst2(s0,sigma(p),Nfft3,gamma);
   
    Nfft4 = 8*N;
    [STFT3,SST3,VSST3_s] = sst2(s0,sigma(p),Nfft4,gamma);
    
    fs = 0:Nfft1/2-1;
    Ren_bis(q,1,p) = renyi_entropy(abs(VSST_s(1:Nfft1/2,:)),t,fs',3);
    fs = 0:Nfft2/2-1; 
    Ren_bis(q,2,p) = renyi_entropy(abs(VSST1_s(1:Nfft2/2,:)),t,fs',3);  
    fs = 0:Nfft3/2-1;
    Ren_bis(q,3,p) = renyi_entropy(abs(VSST2_s(1:Nfft2/2,:)),t,fs',3);
    fs = 0:Nfft4/2-1;
    Ren_bis(q,4,p) = renyi_entropy(abs(VSST3_s(1:Nfft2/2,:)),t,fs',3);
   end 
  end
  Ren(:,:) = mean(Ren_bis);
 end
 figure()
 plot(sigma,Ren(1,:),sigma,Ren(2,:),':',sigma,Ren(3,:),'--',sigma,Ren(4,:),'-.');
end