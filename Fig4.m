
%% Signal
N = 1024;
t = (0:N-1)/N;

a  = 2;
s1 = a.*exp(2*pi*1i*(250*t+50*t.^3));
s2 = a.*exp(2*pi*1i*(130*t+100*t.^2));
s3 = a.*exp(2*pi*1i*(90*t+0.2*cos(3*pi*t)));
s  = s1+s2+s3;
 
%% Parameters
sigma = 0.04;
gamma = 0.01;

% Choice of time and frequency bins
ft = 1:N/2;
bt = 1:N;

% reconstruction parameters
lambda = 0;
beta   = 0;


%% Compute TF transforms
s = add_awgn_noise(s,0); 

%% Compute TF transforms signal s1
%[STFT,SST,VSST]    = sst2_new_gamma_adapt(s,1/sigma^2/N,N);
[STFT,SST,VSST] = sst2(s,sigma,N,gamma);

%Figure A
fs = 0:N/2-1;
imagesc(t,fs,abs(VSST(1:N/2,:)));
set(gca,'YDir','normal')
colormap(1-gray)

Nfft1 = 8*N;
%[STFT,SST,VSST] = sst2_new_gamma_adapt(s,1/sigma^2/N,Nfft1);
[STFT,SST,VSST]  = sst2(s,sigma,Nfft1,gamma);

nr    = 3; %number of modes to be detected
clwin = 10;   
[Cs, Es] = exridge_mult(VSST,nr,lambda,beta,clwin);

for k = nr:-1:1
 %numerical integration of R_vsst
 integ = cumtrapz(t,N/Nfft1*(Cs(k,:)-1));
 sp = s.*exp(-2*1i*pi*(integ-100.*t));
 
 %perform the sunchrosqueezing transform on the demodulated modes 
 %[STFT1,SST1,VSST1] = sst2_new_gamma_adapt(sp,1/sigma^2/N,N);
  [STFT1,SST1,VSST1] = sst2(sp,sigma,N,gamma);

  %Figure B C D
  figure()
  imagesc(t,fs,abs(VSST1(1:N/2,:)));
  if k == 1, 
   set(gca,'ylim',[0 200]);
  elseif k == 2,
    set(gca,'ylim',[0 250]);
  end  
  set(gca,'YDir','normal')
  colormap(1-gray)
end

