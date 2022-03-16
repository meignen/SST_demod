function [sp1_s1,integ1] = demod(s1,VSST_s1,Nfft1,t,N,jump)   
%Computation of the ridges with no regularization, effect of zero padding 
 beta   = 0;
 lambda = 0;
 %for signal s1
 [Cs2, ~] = exridge(VSST_s1,lambda,beta,jump);
 
 R_vsst0  = Cs2-1;
 
 %numerical integration of R_vsst0 and demodulated signal
 integ1 = cumtrapz(t,N/Nfft1*R_vsst0);
 sp1_s1 = s1.*exp(-2*1i*pi*(integ1-100.*t));
end 