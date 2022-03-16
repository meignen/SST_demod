function [STFT,SST,VSST,omega,omega2] = sst2(s,sigma,Nfft,gamma)

 %% sst2_new : computes the STFT of a signal and different versions of synchrosqueezing
 %
 % INPUTS:   
 %   s: real or complex signal
 %   sigma: the variance of the Gaussian window
 %   Nfft: number of frequency bins
 %   gamma: threshold on the STFT for reassignment
 % OUTPUTS:   
 %   STFT : the short-time Fourier transform
 %   SST  : standard synchrosqueezing
 %   VSST1: vertical second-order synchrosqueezing [1]
 % REFERENCES
 % [1] Behera, R., Meignen, S., & Oberlin, T. (2015). Theoretical Analysis
 % of the Second-order Synchrosqueezing Transform. To appear in ACHA
 
 s = s(:);
 N = length(s);          
 
 ft   = 1:Nfft;
 bt   = 1:N;
  
 prec = 10^(-3);
 L =  sigma*N;
 l = floor(L*sqrt(-log(prec)/pi))+1;
 g = amgauss(2*l+1,l+1,L);
 
 % Window definition
  
  n   = (0:2*l)'-l;
  t0  = n/N;
  t0  = t0(:);
  a   = pi/sigma^2;
  gp  = -2*a*t0.*g; 
  gpp = (-2*a+4*a^2*t0.^2).*g; % g''

 % Initialization
 STFT  = zeros(Nfft,N);
 SST   = zeros(Nfft,N);
 VSST  = zeros(Nfft,N);
 
 omega  = zeros(Nfft,N);
 tau    = zeros(Nfft,N);
 omega2 = zeros(Nfft,N);
 phipp  = zeros(Nfft,N);
             
 %% Computes STFT and reassignment operators

 for b=1:N
 	% STFT, window g  
 	time_inst = -min([l,b-1]):min([l,N-b]);
    tmp = fft(s(bt(b)+time_inst).*g(l+time_inst+1),Nfft);
 	vg  = tmp(ft);
     
 	% STFT, window xg           
 	tmp = fft(s(bt(b)+time_inst).*(time_inst)'/N.*g(l+time_inst+1),Nfft);
 	vxg = tmp(ft);
       
    % operator Lx (dtau)
	
    tau(:,b) = vxg./vg;
 	
    % STFT, window gp
 	tmp = fft(s(bt(b)+time_inst).*gp(l+time_inst+1),Nfft);
 	vgp = tmp(ft);
    
    % operator omega
    omega(:,b) = N/Nfft*(ft-1)'-real(vgp/2/1i/pi./vg);    
 	
    
    % STFT, window gpp
 	tmp  = fft(s(bt(b)+time_inst).*gpp(l+time_inst+1),Nfft);
 	vgpp = tmp(ft);
       
    %STFT, windox xgp
 	tmp  = fft(s(bt(b)+time_inst).*(time_inst)'/N.*gp(l+time_inst+1),Nfft);
 	vxgp = tmp(ft);
    
       
 	%computation of the two different omega 
        
    phipp(:,b) = 1/2/1i/pi*(vgpp.*vg-vgp.^2)./(vxg.*vgp-vxgp.*vg);
       
    %new omega2
    omega2(:,b) = omega(:,b) - real(phipp(:,b)).*real(tau(:,b))...
                              + imag(phipp(:,b)).*imag(tau(:,b)); 

	% Storing STFT       
    STFT(:,b) = vg.*exp(2*1i*pi*(ft-1)'*min(l,b-1)/Nfft);%renormalized so that it fits with recmodes
 end
  
 %% reassignment step
 for b=1:N
    for eta=1:Nfft
        if abs(STFT(eta,b))> gamma
           k = 1+round(Nfft/N*omega(eta,b));
            if (k >= 1) && (k <= Nfft)
             % original reassignment
             SST(k,b) = SST(k,b) + STFT(eta,b);
            end
            %reassignment using new omega2
            k = 1+round(Nfft/N*omega2(eta,b));
            if k>=1 && k<=Nfft
                % second-order Vertical reassignment: VSST
                VSST(k,b) = VSST(k,b) + STFT(eta,b);
            end 
        end
    end
 end
end