clc;
clear all;
close all;
[signal fs] = audioread('sp11(1).wav');
 Signal = signal(1:22500)';
nse = audioread('train.wav');
 noise = nse(1:22500);
num = numel(Signal);      %quantity of input signal samples
num1 = numel(noise);      %quantity of noise samples
m = ceil(num/num1);
snr = 15;
noi = repmat(noise,m,1);  %array of noise with length not less than num
noi = noi(1:num);         %array of noise with length equals to num
sa = sqrt(sum(Signal.^2));%numerator of power ratio
sb1 = sqrt(sum(noi.^2));  %denominator of power ratio
varargout(1) = {20*log10(sa/sb1)}; %SNR of SIGNAL corrupted by NOI
ns = sa*10^(-snr/20)/sb1; %coefficient intended to get a specified SNR
s = noi*ns+(Signal)';        %array of corrupted signal with a specified SNR
s=Signal'+noise;
sig=(s);
 
 wlen = 512; timestep = wlen/2/2; numfreq = 320;
 awin = hamming(wlen ) ;
 
Y1 = s(1:16000);
Fo=8000;
nsamp1 = length(Y1); 
numtime1 = ceil((nsamp1-wlen+1)/timestep )-1 ;
tf111 = zeros(numfreq , numtime1);
for i = 1:numtime1
sind1 = ((i-1)*timestep)+1;
tf111(:,i) = fft(Y1(sind1:(sind1+wlen-1)).*awin,numfreq );
end
tf111;
tf11 = tf111(1:(numfreq/2),1:numtime1); 
Signal=Signal';
Fo = 8000;
nsamp2 = length(Y1); 
numtime2 = ceil((nsamp2-wlen+1)/timestep )-1 ;
tf222 = zeros(numfreq , numtime2);
for i2 = 1:numtime2
sind2 = ((i2-1)*timestep)+1;
tf222(:,i2) = fft(Signal(sind2:(sind2+wlen-1)).*awin,numfreq );
end
tf222;
tf222 = tf222(1:(numfreq/2),1:numtime1);
Fo = 8000;
nsamp21 = length(Y1); 
numtime21 = ceil((nsamp21-wlen+1)/timestep )-1 ;
tf121 = zeros(numfreq , numtime21);
for i21 = 1:numtime21
sind21 = ((i21-1)*timestep)+1;
tf121(:,i21) = fft(noise(sind21:(sind21+wlen-1)).*awin,numfreq );
end
tf121;
% tf121=tf121(1:(numfreq/2),1:numtime1);
 
 
swin=sqrt(2)*awin /wlen; 
 swin = swin ( : ) ; % make synthesis window go columnwise
 winlen = length ( swin ) ;
 [numfreq numtime] = size(tf111);
 ind3 = rem((1 : winlen )-1, numfreq)+1;
xN3 = zeros ((numtime-1)* timestep + winlen , 1 ) ;
 for I3 = 1:numtime-1 % overlap , window, and add
 temp3 = (numfreq)*real(ifft(tf121(:,I3)));
 sind3 = ((I3-1)*timestep ) ;
 rind3 = ( sind3 +1):( sind3+winlen ) ;
 xN3(rind3 ) = xN3( rind3 ) + temp3( ind3 ).*swin ;
  %Est(:,I) = xN(1: length (X11));
 end
 
TF=abs(tf11);
 
numtime2 = numtime1;
 
YTT = reshape(TF,(numfreq/2)*numtime1,1);
V = TF;
%V=(abs(TF).^2);
Y = V;
VN = [];
VV1 = sum(V);
for LK = 1:numtime2
VN(:,LK) = V(:,LK)./VV1(:,LK);
end
Vk = VN;
V = Vk./sum(VV1);
surf(log(abs(tf11)));view(0,90);shading interp
 
A111 = rand(numfreq/2,100);
X111 = rand(100,numtime2);
A222 = rand(numfreq/2,100);
X222 = rand(100,numtime2);
A11 = sum(A111);
A22 = sum(A222);
X11 = sum(X111,2);
X22 = sum(X222,2);
A1 = [];
A2 = [];
 
for U = 1:100
    A1 = [A1 A111(:,U)./A11(:,U)];
    
    A2 = [A2 A222(:,U)./A22(:,U)];
    
end
 
X1 = [];
X2 = [];
for UU=1:100
    X1 = [X1; X111(UU,:)./X11(UU,:)];
    X2 = [X2;X222(UU,:)./X22(UU,:)];
end
 
 
% 
A = [A1 A2];
X = [X1;X2];
 
LMD111 = (abs(tf222./std(s))).^0.5; 
 
LMD1 = exp(-LMD111);
    
  LMD222 = (abs(tf121./std(s))).^0.5;
LMD2 = exp(-LMD222);
LMD2 = LMD2(1:160,:)
for k2 = 1:5
    V1 = V.*LMD1;
    V2 = V.*LMD2;
    
   Gam=((A1*X1).*LMD1)+((A2*X2).*LMD2);
   
   M1= V1./Gam;
   M2=V2./Gam;
   A1= A1.*((M1*X1.')./(ones(numfreq/2,numtime2)*X1.'));
   
   A2= A2.*((M2*X2.')./(ones(numfreq/2,numtime2)*X2.'));
   X1=X1.*(A1.'*M1);
   X2=X2.*(A2.'*M2);
end
% 
for  k3=1:5
    F1=(A1*X1)./(A*X);
    
    F2=(A2*X2)./(A*X);
    VCAP1=F1.*V;
        VCAP2=F2.*V;
end
 
E1=V-VCAP1;
E2=V-VCAP2;
 
TFF11=VCAP1(1:numfreq/2,1:numtime1);
 
TF1neg=conj(TFF11(2:(numfreq/2)-1,:));
 
TF1negflp=flipdim(TF1neg,1);
TFF1=[TFF11;TF1negflp];
 
TFF22=VCAP2(1:(numfreq/2),1:numtime1);
 
TF2neg=conj(TFF22(2:(numfreq/2)-1,:));
 
TF2negflp=flipdim(TF2neg,1);
TFF2=[TFF22;TF2negflp];
 
 
 
[DF DT]=size(TFF1);
 ma=angle(tf111);
 ma=ma(1:318,:);
S2=TFF2.*(cos(ma)+sqrt(-1)*(sin(ma)));
 
 
[DF DT]=size(TFF1);
 ma1=angle(tf111);
 ma1=ma1(1:318,:);
S33=TFF1.*(cos(ma1)+sqrt(-1)*(sin(ma1)));
 
 
nsamp1 = length(Y1); 
numtime1 = ceil((nsamp1-wlen+1)/timestep )-1 ;
tf111 = zeros(numfreq , numtime1);
for i = 1:numtime1
sind1 = ((i-1)*timestep)+1;
tf111(:,i) = fft(noise(sind1:(sind1+wlen-1)).*awin,numfreq );
end
tf121;
 
 
 
 
 swin=sqrt(2)*awin /wlen; 
 swin = swin ( : ) ; % make synthesis window go columnwise
 winlen = length ( swin ) ;
 [numfreq numtime] = size(tf11);
ind5 = rem((1 : winlen )-1, numfreq)+1;
xN5 = zeros ((numtime-1)* timestep + winlen , 1 ) ;
 for I5 = 1:numtime-1 % overlap , window, and add
 temp5 = numfreq*real(ifft(tf121(:,I5)));
 sind5 = ((I5-1)*timestep ) ;
 rind5 = ( sind5 +1):( sind5+winlen ) ;
 xN5(rind5 ) = xN5( rind5 ) + temp5( ind5 ).*swin ;
  %Est(:,I) = xN(1: length (X11));
 end
 
 
 swin=sqrt(2)*awin /wlen; 
 swin = swin ( : ) ; % make synthesis window go columnwise
 winlen = length ( swin ) ;
 [numfreq numtime] = size(TFF1);
 ind6 = rem((1 : winlen )-1, numfreq)+1;
xN6 = zeros ((numtime-1)* timestep + winlen , 1 ) ;
 for I6 = 1:numtime-1 % overlap , window, and add
 temp6 = numfreq*real(ifft(S2(:,I6)));
 sind6 = ((I6-1)*timestep ) ;
 rind6 = ( sind6 +1):( sind6+winlen ) ;
 xN6(rind6 ) = xN6( rind6 ) + temp6( ind6 ).*swin ;
  %Est(:,I) = xN(1: length (X11));
 end
 
 for TT=1:12000
 Mse2=10*log10(sqrt(sum(Signal(TT).^2)./12000)./(sqrt(sum((Signal(TT)-xN6(TT))).^2)./12000)); %Mse = Mean square Error
%out(TT)=10*log10(Mse(TT));
 end
 
subplot(3,1,1)
surf(log(abs(tf111)));view(0,90);shading interp
xlabel('Noisy speech');
subplot(3,1,2)
surf(log(abs(S33)));view(0,90);shading interp
subplot(3,1,3)
surf(log(abs(S2)));view(0,90);shading interp


