function [MSE_stft,MSE_sst,MSE_vsst] = Fig2(cas,Nois,SNR)
 %if cas = 1: signal1, if cas = 2: signal2
 %if Nois = 1: noise 
 
 N = 1024;
 t = (0:N-1)/N;

 if ( cas == 1)
  %test signal 1
  a = 2;
  s = a.*exp(2*pi*1i*(130*t+100*t.^2));
  phi1 = 130+200*t;
 else
  %test signal 2   
  a = 2;
  s = a.* exp(2*pi*1i*(330*t+16*cos(3*pi*t)));
  phi1 = 330-48*pi*sin(3*pi*t);
 end
 
 if (Nois == 0)
  %% Parameters
  gamma = 0.01;
  sigma = 0.04;

  %% Compute TF transforms
  
  Nfft1 = 8*N;
 [STFT,SST,VSST] = sst2(s,sigma,Nfft1,gamma); 
 
  %classical ridge extraction
  lambda = 0:0.0001:0.002;
  beta   = 0:0.0002:0.004;
  
  MSE_stft = zeros(length(lambda),length(beta));
  MSE_sst = zeros(length(lambda),length(beta));
  MSE_vsst = zeros(length(lambda),length(beta));
  index = 100:N-100;
  
  for k = 1:length(lambda)
   for p = 1:length(beta)   
   
    [Cs0,~] = exridge(STFT,lambda(k),beta(p),10); 
    
    [Cs1,~] = exridge(SST,lambda(k),beta(p),10);
  
    [Cs2,~] = exridge(VSST,lambda(k),beta(p),10);
    
    MSE_stft(k,p) =  sqrt(sum((N/Nfft1*(Cs0(index)-1)-phi1(index)).^2)/(length(index)-1));
    MSE_sst(k,p)  =  sqrt(sum((N/Nfft1*(Cs1(index)-1)-phi1(index)).^2)/(length(index)-1));
    MSE_vsst(k,p) =  sqrt(sum((N/Nfft1*(Cs2(index)-1)-phi1(index)).^2)/(length(index)-1));
    end
  end 
 end
 
 if (Nois == 1)
  lambda = 0:0.0001:0.002;
  beta   = 0:0.0002:0.004;
  MSE_stft_iter = zeros(30,length(lambda),length(beta));
  MSE_sst_iter  = zeros(30,length(lambda),length(beta));
  MSE_vsst_iter = zeros(30,length(lambda),length(beta));
  %% Parameters
  gamma = 0.01;
  sigma = 0.04;
  index = 100:N-100;
 
  %% Compute TF transforms
  
  %classical ridge extraction
  for q = 1:30
   s0 = add_awgn_noise(s,SNR);
   Nfft1 = 8*N;
   [STFT,SST,VSST] = sst2(s0,sigma,Nfft1,gamma); 
   
 
   for k = 1:length(lambda)
    for p = 1:length(beta)   
   
     [Cs0,~] = exridge(STFT,lambda(k),beta(p),10); 
     [Cs1,~] = exridge(SST,lambda(k),beta(p),10); 
     [Cs2,~] = exridge(VSST,lambda(k),beta(p),10);
  
     MSE_stft_iter(q,k,p) =  sqrt(sum((N/Nfft1*(Cs0(index)-1)-phi1(index)).^2)/(length(index)-1));
     MSE_sst_iter(q,k,p)  =  sqrt(sum((N/Nfft1*(Cs1(index)-1)-phi1(index)).^2)/(length(index)-1));
     MSE_vsst_iter(q,k,p) =  sqrt(sum((N/Nfft1*(Cs2(index)-1)-phi1(index)).^2)/(length(index)-1));
    end
   end 
  end
  
  MSE_stft = zeros(length(lambda),length(beta));
  MSE_sst  = zeros(length(lambda),length(beta));
  MSE_vsst = zeros(length(lambda),length(beta));
  
  X =  mean(MSE_stft_iter);
  Y =  mean(MSE_sst_iter);
  Z =  mean(MSE_vsst_iter);
  
  MSE_stft(:,:) = X(1,:,:);
  MSE_sst(:,:)  = Y(1,:,:);
  MSE_vsst(:,:) = Z(1,:,:);
 end