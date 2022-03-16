function [sp1_s,sp2_s,sp3_s,sp4_s,integ1,integ2,integ3,integ4] = ...
 demod_multi(s,VSST_s,VSST1_s,VSST2_s,VSST3_s,Nfft1,Nfft2,Nfft3,Nfft4,t,N,nr,jump)   

%Computation of the ridges with no regularization, effect of zero padding 
 beta   = 0;
 lambda = 0;
 
 %computation of the different ridges
 [Cs2, ~] = exridge_mult(VSST_s,nr,lambda,beta,jump);
 [Cs21,~] = exridge_mult(VSST1_s,nr,lambda,beta,jump);
 [Cs22,~] = exridge_mult(VSST2_s,nr,lambda,beta,jump);
 [Cs23,~] = exridge_mult(VSST3_s,nr,lambda,beta,jump);
 
 integ1 = zeros(size(Cs2));
 integ2 = zeros(size(Cs21));
 integ3 = zeros(size(Cs22));
 integ4 = zeros(size(Cs23));
 sp1_s  =  zeros(nr,length(s));
 sp2_s  =  zeros(nr,length(s));
 sp3_s  =  zeros(nr,length(s));
 sp4_s  =  zeros(nr,length(s));
 
 for k = 1:nr
  %numerical integration of R_vsst0 and demodulated signal
  integ1(k,:) = cumtrapz(t,N/Nfft1*(Cs2(k,:)-1));
  sp1_s(k,:)  = s.*exp(-2*1i*pi*(integ1(k,:)-100.*t));
 
  %numerical integration of R_vsst1
  integ2(k,:) = cumtrapz(t,N/Nfft2*(Cs21(k,:)-1));
  sp2_s(k,:)  = s.*exp(-2*1i*pi*(integ2(k,:)-100.*t));
 
  %numerical integration of R_vsst2
  integ3(k,:) = cumtrapz(t,N/Nfft3*(Cs22(k,:)-1));
  sp3_s(k,:)  = s.*exp(-2*1i*pi*(integ3(k,:)-100.*t));
  
  %numerical integration of R_vsst3
  integ4(k,:) = cumtrapz(t,N/Nfft4*(Cs23(k,:)-1));
  sp4_s(k,:)  = s.*exp(-2*1i*pi*(integ4(k,:)-100.*t));
  
 end 