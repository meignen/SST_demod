function [snr_direct,snr_demod] = Fig7_8(cas)
 
 N = 1024;

 % Parameters
 gamma = 10^(-3);
 index = N/8+1:7*N/8;
 lambda = 0;
 beta   = 0;

 t = (0:N-1)/N;
 sigma = 0.01:0.01:0.2;
 
 %% Tests signals 1
 if (cas == 1)
  a  = 2;
  s1 = a.*exp(2*pi*1i*(250*t+50*t.^3));
  s2 = a.*exp(2*pi*1i*(130*t+100*t.^2));
  s  = s1+s2;
  clwin = 10;
  nr = 2;
 else
 %% Test signal 2
  a  = 2;
  s1 = a.*exp(2*pi*1i*(330*t+16*cos(3*pi*t)));
  s2 = a.*exp(2*pi*1i*(190*t+9*cos(3*pi*t)));
  s  = s1+s2;
  clwin = 30;
  nr = 2;
 end

 %% computation of the optimal sigma
  fs = 0:N/2-1;
  for p = 1:length(sigma)
   [~,~,VSST_s] = sst2(s,sigma(p),N/2,gamma);
   E = renyi_entropy(abs(VSST_s),t,fs',3);
   if p == 1
    En = E;
    p_opt = p;
   else
    if (E < En)
     En = E;   
     p_opt = p;
    end
   end  
  end

 sigma_opt = sigma(p_opt);
 
 %%computation of transforms with different frequency resolutions 
 Nfft1 = N;
 [STFT,SST,VSST_s]    = sst2(s,sigma_opt,Nfft1,gamma);

 Nfft2 = 2*N;
[STFT1,SST1,VSST1_s] = sst2(s,sigma_opt,Nfft2,gamma);

 Nfft3 = 4*N;
[STFT2,SST2,VSST2_s] = sst2(s,sigma_opt,Nfft3,gamma);

 Nfft4 = 8*N;
[STFT3,SST3,VSST3_s] = sst2(s,sigma_opt,Nfft4,gamma);

 [sp1_s,sp2_s,sp3_s,sp4_s,integ1,integ2,integ3,integ4] = ...
 demod_multi(s,VSST_s,VSST1_s,VSST2_s,VSST3_s,Nfft1,Nfft2,Nfft3,Nfft4,t,N,nr,clwin);
 
%ground truth
 if (cas == 1)
  integ1_true    = cumtrapz(t,250+150*t.^2);
  integ2_true    = cumtrapz(t,130+200*t);
  gt1            = s.*exp(-2*1i*pi*(integ1_true-100.*t));
  gt2            = s.*exp(-2*1i*pi*(integ2_true-100.*t));
  gt= [gt1;gt2];
  integ_true = [integ1_true;integ2_true]; 
 else
  integ1_true    = cumtrapz(t,330-16*3*pi*sin(3*pi*t));
  integ2_true    = cumtrapz(t,190-9*3*pi*sin(3*pi*t));
  gt1            = s.*exp(-2*1i*pi*(integ1_true-100.*t));
  gt2            = s.*exp(-2*1i*pi*(integ2_true-100.*t));
  gt= [gt1;gt2];
  integ_true = [integ1_true;integ2_true];  
 end
 
 %%computation of the ridges for direct reconstruction
 [Cs2,~] = exridge_mult(VSST_s,nr,lambda,beta,clwin);

 %%signal reconstruction
 d = 0:1:10;
 snr_direct = zeros(nr,length(d));
 snr_demod_true = zeros(nr,length(d));
 snr_demod  = zeros(4,nr,length(d)); %demodulated signal depending on frequency resolution

 signal = [s1;s2];
 sign1  = zeros(nr,length(s));
 
 for  k = 1:length(d) %influence of parameter d
  k
  for p = 1:nr
   sign1(p,:) = 1/Nfft1*recmodes(VSST_s,Cs2(p,:),d(k));   
   
   snr_direct(p,k) =  snr(signal(p,index),sign1(p,index)-signal(p,index));
   
   [~,~,VSST_sd_true] = sst2(gt(p,:),sigma_opt,Nfft1,gamma);
   [C,~]                = exridge_mult(VSST_sd_true,p,lambda,beta,clwin);  
   imf                  = 1/Nfft1*recmodes(VSST_sd_true,C(p,:),d(k));
   imf                  = imf.*exp(2*1i*pi*(integ_true(p,:)-100.*t));
   snr_demod_true(p,k)  = snr(signal(p,index),imf(index)-signal(p,index));
   
   [~,~,VSST_sd_1] = sst2(sp1_s(p,:),sigma_opt,Nfft1,gamma);
   [C,~]                  = exridge_mult(VSST_sd_1,p,lambda,beta,clwin);  
   imf                    = 1/Nfft1*recmodes(VSST_sd_1,C(p,:),d(k));
   imf                    = imf.*exp(2*1i*pi*(integ1(p,:)-100.*t));
   snr_demod(1,p,k)       = snr(signal(p,index),imf(index)-signal(p,index));
   
   [~,~,VSST_sd_1] = sst2(sp2_s(p,:),sigma_opt,Nfft1,gamma);
   [C,~]                  = exridge_mult(VSST_sd_1,p,lambda,beta,clwin);  
   imf                    = 1/Nfft1*recmodes(VSST_sd_1,C(p,:),d(k));
   imf                    = imf.*exp(2*1i*pi*(integ2(p,:)-100.*t));
   snr_demod(2,p,k)       = snr(signal(p,index),imf(index)-signal(p,index));
   
   [~,~,VSST_sd_1] = sst2(sp3_s(p,:),sigma_opt,Nfft1,gamma);
   [C,~]                  = exridge_mult(VSST_sd_1,p,lambda,beta,clwin);  
   imf                    = 1/Nfft1*recmodes(VSST_sd_1,C(p,:),d(k));
   imf                    = imf.*exp(2*1i*pi*(integ3(p,:)-100.*t));
   snr_demod(3,p,k)       = snr(signal(p,index),imf(index)-signal(p,index));
   
   [~,~,VSST_sd_1] = sst2(sp4_s(p,:),sigma_opt,Nfft1,gamma);
   [C,~]                = exridge_mult(VSST_sd_1,p,lambda,beta,clwin);  
   imf                    = 1/Nfft1*recmodes(VSST_sd_1,C(p,:),d(k));
   imf                    = imf.*exp(2*1i*pi*(integ4(p,:)-100.*t));
   snr_demod(4,p,k)       = snr(signal(p,index),imf(index)-signal(p,index));

  end
 end
 
 figure()
 fs = N/Nfft1*(0:Nfft1/2-1);
 imagesc(t,fs,abs(VSST_s(1:Nfft1/2,:)));
 colormap(1-gray)
 set(gca,'ydir','normal');
 
 Y=zeros(1,length(d));
 Z=zeros(1,length(d));

 figure()
 Y(:)  = snr_demod(1,2,:);
 Z(:)  = snr_demod(1,1,:);
 plot(d,snr_demod_true(2,:),d,snr_demod_true(1,:),'--',... 
      d,Y,'-*',d,Z,'-d');%,...
 
 figure()
 Y(:) = snr_demod(3,2,:);
 Z(:) = snr_demod(3,1,:);
 plot(d,snr_demod_true(2,:),d,snr_demod_true(1,:),':',... 
      d,Y,'-*',d,Z,'-d');%,...

 figure()
 Y(:) = snr_demod(4,2,:);
 Z(:) = snr_demod(4,1,:);
 plot(d,snr_demod_true(2,:),d,snr_demod_true(1,:),'--',... 
      d,snr_direct(2,:),':',d,snr_direct(1,:),'-.',... 
      d,Y,'-*',d,Z,'-d');%,...

end
 