 
%% Signal
N = 1024;
t = (0:N-1)/N;

a  = 2;
s1 = a.*exp(2*pi*1i*(250*t+50*t.^3));
s2 = a.*exp(2*pi*1i*(130*t+100*t.^2));
s3 = a.*exp(2*pi*1i*(90*t+0.2*cos(3*pi*t)));
jump = 10;

%% Parameters
 gamma = 0.01;
 sigma = 0.04;

 % Choice of time and frequency bins
 ft = 1:N/2;
 bt = 1:N;
 
 %% Compute TF transforms signal s1
 Nfft1 = N;
[~,~,VSST_s1]    = sst2(s1,sigma,Nfft1,gamma);
 
 %% Compute TF transforms signal s2
 Nfft1 = N;
[~,~,VSST_s2]     = sst2(s2,sigma,Nfft1,gamma);
 
 %% Compute TF transforms signal s3
 Nfft1 = N;
[~,~,VSST_s3]     = sst2(s3,sigma,Nfft1,gamma);
 
 %the VSST Figure 3 A B C 
 t  = (0:N-1)/N;
 fs = N/Nfft1*(0:Nfft1/2-1);
 imagesc(t,fs,abs(VSST_s2(1:Nfft1/2,:)));set(gca,'ydir','normal');set(gca,'ylim',[100 350]);colormap(1-gray);
 figure()
 imagesc(t,fs,abs(VSST_s1(1:Nfft1/2,:)));set(gca,'ydir','normal');set(gca,'ylim',[200 512]);colormap(1-gray);
 figure()
 imagesc(t,fs,abs(VSST_s3(1:Nfft1/2,:)));set(gca,'ydir','normal');set(gca,'ylim',[80 100]);colormap(1-gray);
 
 %with noise
 s1 = add_awgn_noise(s1,0);
 s2 = add_awgn_noise(s2,0);
 s3 = add_awgn_noise(s3,0);
 
 %% Compute TF transforms signal s1
 Nfft1 = N;
[STFT,SST,VSST_s1]    = sst2(s1,sigma,Nfft1,gamma);
 
 %% Compute TF transforms signal s2
[STFT,SST,VSST_s2]     = sst2(s2,sigma,Nfft1,gamma);

 %% Compute TF transforms signal s3
[STFT,SST,VSST_s3]     = sst2(s3,sigma,Nfft1,gamma);

 %the VSST Figure 3 D  E F
 t  = (0:N-1)/N;
 fs = N/Nfft1*(0:Nfft1/2-1);
 figure()
 imagesc(t,fs,abs(VSST_s2(1:Nfft1/2,:)));set(gca,'ydir','normal');set(gca,'ylim',[100 350]);colormap(1-gray);
 figure()
 imagesc(t,fs,abs(VSST_s1(1:Nfft1/2,:)));set(gca,'ydir','normal');set(gca,'ylim',[200 512]);colormap(1-gray);
 figure()
 imagesc(t,fs,abs(VSST_s3(1:Nfft1/2,:)));set(gca,'ydir','normal');set(gca,'ylim',[80 100]);colormap(1-gray);
 
 %computation of demodulated signals 
 %% Compute TF transforms signal s1
 Nfft1 = 8*N;
[STFT,SST,VSST_s1]    = sst2(s1,sigma,Nfft1,gamma);
 
 %% Compute TF transforms signal s2
[STFT,SST,VSST_s2]     = sst2(s2,sigma,Nfft1,gamma);

 %% Compute TF transforms signal s3
[STFT,SST,VSST_s3]     = sst2(s3,sigma,Nfft1,gamma);

%  [sp1_s1,sp2_s1,sp3_s1,sp4_s1,sp5_s1] = demod(s1,VSST_s1,VSST1_s1,VSST2_s1,VSST3_s1,VSST4_s1,...
%                                                       Nfft1,Nfft2,Nfft3,Nfft4,Nfft5,t,N,jump);   
% 
%  [sp1_s2,sp2_s2,sp3_s2,sp4_s2,sp5_s2] = demod(s2,VSST_s2,VSST1_s2,VSST2_s2,VSST3_s2,VSST4_s2,...
%                                                       Nfft1,Nfft2,Nfft3,Nfft4,Nfft5,t,N,jump);
%  
%  [sp1_s3,sp2_s3,sp3_s3,sp4_s3,sp5_s3] = demod(s3,VSST_s3,VSST1_s3,VSST2_s3,VSST3_s3,VSST4_s3,...
%                                                       Nfft1,Nfft2,Nfft3,Nfft4,Nfft5,t,N,jump);                                                

 [sp1_s1] = demod(s1,VSST_s1,Nfft1,t,N,jump);   
 [sp1_s2] = demod(s2,VSST_s2,Nfft1,t,N,jump);
 [sp1_s3] = demod(s3,VSST_s3,Nfft1,t,N,jump);                                                

%computation of the VSST of the demodulated modes
 %Figure 3 G H I
 figure()
 [STFT,SST,VSST]    = sst2(sp1_s1,sigma,N,gamma);
 imagesc(t,fs,abs(VSST(1:N/2,:)));
 set(gca,'ydir','normal');set(gca,'ylim',[80 120]);colormap(1-gray);
 
 figure()
 [STFT,SST,VSST]    = sst2(sp1_s2,sigma,N,gamma);
 imagesc(t,fs,abs(VSST(1:N/2,:)));
 set(gca,'ydir','normal');set(gca,'ylim',[80 120]);colormap(1-gray);
 
 figure()
 [STFT,SST,VSST]    = sst2(sp1_s3,sigma,N,gamma);
 imagesc(t,fs,abs(VSST(1:N/2,:)));
 set(gca,'ydir','normal');
 set(gca,'ylim',[80 120]);
 colormap(1-gray);
 
 
         
 
 

