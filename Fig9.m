function [snr_direct,snr_demod] = Fig10(cas)
 
 N = 1024;

 % Parameters
 index = N/8+1:7*N/8;
 lambda = 0;
 beta   = 0;

 t = (0:N-1)/N;
 %sigma = 0.01:0.01:0.2;
 gamma = 10^(-3);
 
 %% Tests signals 1
 if (cas == 1)
  a  = 2;%exp(-8*(t-0.5).^2);
  s1 = a.*exp(2*pi*1i*(250*t+50*t.^3));
  s2 = a.*exp(2*pi*1i*(130*t+100*t.^2));
  %s3 = a.*exp(2*pi*1i*(90*t+0.2*cos(3*pi*t)));
  s  = s1+s2;%+s3;
  nr = 2;
  clwin = 10;
  SNR = -5:15;
 else
 %% Test signal 2
  a  = 2;%exp(-8*(t-0.5).^2);
  s1 = a.*exp(2*pi*1i*(330*t+16*cos(3*pi*t)));
  s2 = a.*exp(2*pi*1i*(190*t+9*cos(3*pi*t)));
  %s3 = a.*exp(2*pi*1i*(40*t));
  s  = s1+s2;%+s3;
  nr = 2;
  clwin = 30;
  SNR = 0:15;
 end
 
 d   = 0:5:5;
 
 signal = [s1;s2];
 snr_direct = zeros(nr,length(SNR),length(d));
 snr_demod  = zeros(nr,length(SNR),length(d));

 snr_direct_int = zeros(10,nr,length(SNR),length(d));
 snr_demod_int  = zeros(10,nr,length(SNR),length(d));
 
 for r = 1:10
  r
  for q = 1:length(SNR)
   q
   s0 = add_awgn_noise(s,SNR(q));  
   %computation of the optimal sigma
%    if (SNR(q) >= 0)
%       fs = 0:N/2-1;
%       for p = 1:length(sigma)
%       [STFT,SST,VSST_s] = sst2_new_gamma_adapt(s0,1/sigma(p)^2/N,N);
%       E = renyi_entropy(abs(VSST_s),t,fs',3);
%       if p == 1, 
%        En = E;
%        p_opt = p;
%       else
%        if (E < En)
%         En = E;   
%         p_opt = p;
%        end
%       end  
%      end
%     else
%      En = 0;
%       for p = 1:length(sigma)
%      [STFT,SST,VSST_s] = sst2_new_gamma_adapt(s0,1/sigma(p)^2/N,N,N);
%      [Cs2, Es] = exridge_mult(VSST_s,nr,lambda,beta,clwin);
%      E = sum(Es)/sum((abs(VSST_s(:))).^2);
%      if (E > En)
%       En = E;
%       p_opt = p;
%      end 
%     end
%    end
   
   %sigma_opt = sigma(p_opt);
   
   sigma_opt = 0.04;
   gamma = 10^(-6);
   %%computation of demodulated signals+instantaneous frequencies
   Nfft1 = 8*N;
   %[STFT3,SST3,VSST_s] = sst2_new_gamma_adapt(s0,1/sigma_opt^2/N,Nfft1);
   [STFT3,SST3,VSST_s] = sst2(s0,sigma_opt,Nfft1,gamma);
   [Cs2,~] = exridge_mult(VSST_s,nr,lambda,beta,clwin);
   integ1 = zeros(size(Cs2));
   sp1_s  =  zeros(nr,length(s));
   for k = 1:nr
    %numerical integration and demodulated signal
    integ1(k,:) = cumtrapz(t,N/Nfft1*(Cs2(k,:)-1));
    sp1_s(k,:)  = s0.*exp(-2*1i*pi*(integ1(k,:)-100.*t));
   end
   
   %%computation of the ridges for direct reconstruction
   %[STFT3,SST3,VSST_s] = sst2_new(s0,1/sigma_opt^2/N,N,gamma);
   [STFT3,SST3,VSST_s] = sst2(s0,sigma_opt,N,gamma);
   [Cs2,~] = exridge_mult(VSST_s,nr,lambda,beta,clwin);
 
   sign1  = zeros(nr,length(s));
 
   for p = 1:nr
    %[STFT1,SST1,VSST_sd] = sst2_new(sp1_s(p,:),1/sigma_opt^2/N,N,gamma);
    [STFT1,SST1,VSST_sd] = sst2(sp1_s(p,:),sigma_opt,N,gamma);
    %we look for the ridge around 100
    VSST_sd_int  = zeros(size(VSST_sd));
    VSST_sd_int(100-clwin:100+clwin,:) = VSST_sd(100-clwin:100+clwin,:);
    [C, ~] = exridge(VSST_sd_int,lambda,beta,clwin);  
    
    for  k = 1:length(d) %influence of parameter d
     sign1(p,:) = 1/N*recmodes(VSST_s,Cs2(p,:),d(k));   
     snr_direct_int(r,p,q,k) =  snr(signal(p,index),sign1(p,index)-signal(p,index));
     
     imf                    = 1/N*recmodes(VSST_sd,C,d(k));
     imf                    = imf.*exp(2*1i*pi*(integ1(p,:)-100.*t));
     snr_demod_int(r,p,q,k) = snr(signal(p,index),imf(index)-signal(p,index));
    end
   end
  end
  snr_direct = mean(snr_direct_int);
  snr_demod  = mean(snr_demod_int); %demodulated signal depending on frequency resolution
 end 
end
