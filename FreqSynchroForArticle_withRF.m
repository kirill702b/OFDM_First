clear all;close all;
FreqOffset = 0.7; % Frequency offcet in subcarrier spacing
SNR = -20;
%
Nc = 1000; % Number of subcarriers
N = 1024; % Number of IFFT points
Rd = 18e6; % Data rate bit per second
CP = 0.1; % Cyclic prefix length as part of symbol period
Fs1 =(Rd/(1-CP));%20mHz, 50ns
% 
NumOfOsc = 16; % Number of Oscillator for Jakes' model
Fc = 2e9; % Hz, Carrier frequency 2GHz
OsFac = 10; %oversampling factor
Fs = Fc*OsFac; %sampling frequency Fs = 20GHz
VehSpeed = 150; % kmph, Speed of vehicle
NumOfPath = 4; % Number of paths
DelOfPath = [0, 1, 2, 3]/Fs; % Delays of the paths in samples
AvRecPwr = [0, -3, -6, -9]; % Average received power, dB
Nnoise = 3000;%number of noise samples
phiRand = 2*pi/3;

% 

frs = (randi(2,[1 N/2])-1.5)*2;%mseq(2,log2(N/2)); последовательность +-1
frs2 = (randi(2,[1 N])-1.5)*2;% вторая последовательность
A = ifft(frs,N/2);%A -последовательность длины N/2
B = ifft(frs2,N);% последовательность длины N
% seq1df = seq1.*
%  seq1 = randn(1,N/2)+1i*randn(1,N/2);
AA = [A(1:N/2),A(1:N/2)];% последовательность длины N
AA_B = [AA, zeros(1,fix(CP*N)),B]; %две последовательности длины 2048+102=2150
% seq1 = [randn(1,3*N)+1i*randn(1,3*N),AA,randn(1,3*N)+1i*randn(1,3*N)];

seq1 = [zeros(1,3*N),real(AA_B),zeros(1,3*N)]+1i* [zeros(1,3*N),imag(AA_B),zeros(1,3*N)];% последовательность с началом и концом без сигнала
% длина 6*N+2150 = 8294

seq1Up = resample(seq1,Fs/Fs1,1); %последовательность на частоте 2ГГц, была на 20MHz
tUp = (1:length(seq1Up))/Fs;% вектор времен для этой новой последовательности 
seq1UpFc = seq1Up.*exp(1i*2*pi*(Fc+FreqOffset*Fs1/N)*tUp);% последовательность уже на несущей + доплеровкий сдвиг

tt1 = (1:length(AA)*Fs/Fs1)/Fs;
nfft = 2^15;
% [pw,f]=pwelch(seq1UpFcAwgn(3*N*Fs/Fs1+1:3*N*Fs/Fs1+length(AA)*Fs/Fs1),[],[],nfft,Fs);
% figure;plot(f,10*log10(pw));grid on;
% 
% [pw,f]=pwelch(seq1UpFcAwgn,[],[],nfft,Fs);
% figure;plot(f,10*log10(pw));grid on;

% aa1 = resample(AA,Fs/Fs1,1).*exp(1i*2*pi*Fc*tt1);
% [pw,f]=pwelch(aa1,[],[],nfft,Fs);
% figure;plot(f-Fs/2,[pw(nfft/2+1:nfft);pw(1:nfft/2)]);grid on;


% seq1 = awgn([zeros(1,3*N),real(AA_B),zeros(1,3*N)],SNR)+1i*awgn( [zeros(1,3*N),imag(AA_B),zeros(1,3*N)],SNR);
%'SampleRate',Fs,
% RayleighCh = comm.RayleighChannel('PathDelays',DelOfPath,'AveragePathGains',AvRecPwr,'FadingTechnique','Sum of sinusoids','NumSinusoids',NumOfOsc);
RayCh1 = rayleighchan(1/Fs,277.3,DelOfPath,AvRecPwr);
% seq1 = step(RayleighCh, (([zeros(1,3*N),real(AA_B),zeros(1,3*N)])+1i*( [zeros(1,3*N),imag(AA_B),zeros(1,3*N)]) ).').';


seq1UpFcAwgn = awgn(seq1UpFc,SNR,'measured');%сигнал после АГБШ канала
% seq1UpFcRay = step(RayleighCh,seq1UpFc.').';%сигнал после Рэлеевского канала
% seq1UpFcRay = filter(RayCh1,seq1UpFc);
seq1UpFcRay = awgn(filter(RayCh1,seq1UpFc),SNR,'measured');
% Ndiv = 10000;
% BufS = fix(length(seq1UpFc)/Ndiv);
% for j =1:Ndiv-1
%     seq1UpFcRay((j-1)*BufS+1:j*BufS) = step(RayleighCh,seq1UpFc((j-1)*BufS+1:j*BufS).').';
% end
% seq1UpFcRay((Ndiv-1)*BufS+1:end) = step(RayleighCh,seq1UpFc((Ndiv-1)*BufS+1:end).').';
seq1UpFcAwgnF0 = seq1UpFcAwgn.*exp(-1i*2*pi*Fc*tUp);%сигнал после переноса частоты на 0
seq1UpFcRayF0 = seq1UpFcRay.*exp(-1i*2*pi*Fc*tUp);%сигнал после переноса частоты на 0
% figure;hold on;grid on;
% plot(seq1UpFcAwgnF0);
% plot(seq1UpFcRayF0,'r');

[pw,f]=pwelch(seq1UpFcAwgnF0,[],[],nfft,Fs);% спектр сигнала 
figure;hold on;grid on;plot(f,10*log10(pw));grid on;%рисуем

[pw,f]=pwelch(seq1UpFcRayF0,[],[],nfft,Fs);% спектр сигнала 
plot(f,10*log10(pw),'r');grid on;%рисуем
%%
seq1UpFcAwgnF0Dec500 = decimate(decimate(decimate(decimate(seq1UpFcAwgnF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
seq1UpFcRayF0Dec500 = decimate(decimate(decimate(decimate(seq1UpFcRayF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)

seq1UpFcAwgnF0Dec500Lpf=seq1UpFcAwgnF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcAwgnF0Dec500);%фильтрация для итоговой децимации
seq1UpFcRayF0Dec500Lpf = seq1UpFcRayF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcRayF0Dec500);% фильтрация для итоговой децимации

seq1UpFcAwgnF0Dec500LpfFs1 = decimate(seq1UpFcAwgnF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты
seq1UpFcRayF0Dec500LpfFs1 = decimate(seq1UpFcRayF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты

for j=1:length(seq1UpFcAwgnF0Dec500LpfFs1)-N %ищем начало (синхронизация) для АГБШ канала
    P1 = sum(conj(seq1UpFcAwgnF0Dec500LpfFs1(1,j:j+N/2-1)).*seq1UpFcAwgnF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2));
    R1 = sum(power(abs(seq1UpFcAwgnF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2)),2));
    RespOfFind(1,j) = power(abs(P1),2)/power(R1,2); % отклик на поиск
end

for j=1:length(seq1UpFcRayF0Dec500LpfFs1)-N %ищем начало (синхронизация) для рэлеевского канала
    P1 = sum(conj(seq1UpFcRayF0Dec500LpfFs1(1,j:j+N/2-1)).*seq1UpFcRayF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2));
    R1 = sum(power(abs(seq1UpFcRayF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2)),2));
    RespOfFind(2,j) = power(abs(P1),2)/power(R1,2); %отклик на поиск
end

tResp = 1:length(RespOfFind);
figure;hold on;grid on;%рисуем отклики во времени
plot(tResp(2.2*N:end-3*N),abs(RespOfFind(1,2.2*N:end-3*N)));%2.2*N:end-3*N
plot(tResp(2.2*N:end-3*N),abs(RespOfFind(2,2.2*N:end-3*N)),'r');
% plot(abs(RespOfFind(1,:)));%2.2*N:end-3*N
% plot(abs(RespOfFind(2,:)),'r');


% for j=1:length(seq1(1,:))-N
%     P1 = sum(conj(seq1(1,j:j+N/2-1)).*seq1(1,j+N/2:j+N/2-1+N/2));
%     R1 = sum(abs(conj(seq1(1,j:j+N/2-1)).*seq1(1,j+N/2:j-1+N/2+N/2)));%sum(abs(conj(seq1(1,j:j+N/2-1)).*seq1(1,j+N/2:j+N/2-1+N/2)));%      %sum(power(abs(seq1(1,j+N/2:j+N/2-1+N/2)),2));
%     RespOfFind(1,j) = power(abs(P1),2)/power(R1,2);
% end

Vk = sqrt(2).*frs2(1:2:end)./frs;%вспомогательная последовательность
% figure;plot(abs(seq1));
t1 = (1:length(AA(1,:)));
t12 = (1:length(AA_B(1,:)));




% seq1df(1,:) = seq1(1,3*N+1:3*N+N).*exp(1i*(2*pi*FreqOffset*t1/N));%Сначала брал только неискаженную последовательность, сдвинутую
% seq12df(1,:) = seq1(1,3*N+1:3*N+N+N+fix(CP*N)).*exp(1i*(2*pi*FreqOffset*t12/N));


seq1df(1,:) = seq1UpFcAwgnF0Dec500LpfFs1(1,3*N+1:3*N+N);%Теперь уже все по честному 
seq12df(1,:) = seq1UpFcAwgnF0Dec500LpfFs1(1,3*N+1:3*N+N+N+fix(CP*N));%АГБШ

seq1df(2,:) = seq1UpFcRayF0Dec500LpfFs1(1,3*N+1:3*N+N);%Теперь уже все по честному
seq12df(2,:) = seq1UpFcRayF0Dec500LpfFs1(1,3*N+1:3*N+N+N+fix(CP*N));%РЭЛЕЙ



FkAwgn = fft(conj(AA).*seq1df(1,:),N);
FkAwgn = [FkAwgn(N/2+1:end),FkAwgn(1:N/2)];
kmax =  find(abs(FkAwgn)==max(abs(FkAwgn)));
FcoarseAwgn = kmax-N/2;
if abs(FkAwgn(kmax-1))<= abs(FkAwgn(kmax+1))
    alp = 1;
else
    alp = -1;
end
FfineAwgn = alp/(abs(FkAwgn(kmax))/abs(FkAwgn(kmax+alp))+1);


FkRay = fft(conj(AA).*seq1df(2,:),N);
FkRay = [FkRay(N/2+1:end),FkRay(1:N/2)];
kmax =  find(abs(FkRay)==max(abs(FkRay)));
FcoarseRay = kmax-N/2;
if abs(FkRay(kmax-1))<= abs(FkRay(kmax+1))
    alp = 1;
else
    alp = -1;
end
FfineRay = alp/(abs(FkRay(kmax))/abs(FkRay(kmax+alp))+1);


% Schmidl
P1 = sum(conj(seq1df(1,1:N/2)).*seq1df(1,1+N/2:N));
phi = angle(P1);
feAwgn = phi/pi;
Corr = exp(-1i*2*phi*t12/N);
seq12dfCorr = seq12df(1,:).*Corr;
F1 = fft(seq12dfCorr(1:N));F1 = [F1(N/2+1:end),F1(1:N/2)];
F2 = fft(seq12dfCorr(N+fix(CP*N)+1:end));F2 = [F2(N/2+1:end),F2(1:N/2)];
F1F2 = conj(F1).*F2;
% B(1) = power(abs(sum(F1F2.*conj(frs(1:2:end)))),2)/2/power(sum(power(F2,2)),2);
for i = 1:N/2
    F1F2sh = circshift(F1F2,[0,2*(i-1)]);
    BgAwgn(i) = power(abs(sum(  F1F2sh(1:2:end).*  conj(Vk)      )),2)/2/power(sum(power(abs(F2),2)),2);
end

P1 = sum(conj(seq1df(2,1:N/2)).*seq1df(2,1+N/2:N));
phi = angle(P1);
feRay = phi/pi;
Corr = exp(-1i*2*phi*t12/N);
seq12dfCorr = seq12df(2,:).*Corr;
F1 = fft(seq12dfCorr(1:N));F1 = [F1(N/2+1:end),F1(1:N/2)];
F2 = fft(seq12dfCorr(N+fix(CP*N)+1:end));F2 = [F2(N/2+1:end),F2(1:N/2)];
F1F2 = conj(F1).*F2;
% B(1) = power(abs(sum(F1F2.*conj(frs(1:2:end)))),2)/2/power(sum(power(F2,2)),2);
for i = 1:N/2
    F1F2sh = circshift(F1F2,[0,2*(i-1)]);
    BgRay(i) = power(abs(sum(  F1F2sh(1:2:end).*  conj(Vk)      )),2)/2/power(sum(power(abs(F2),2)),2);
end

%%
FdProposedAwgn = FcoarseAwgn + FfineAwgn - 1
FdSchmidlAwgn = (1+256-find(BgAwgn==max(BgAwgn)))*2 +feAwgn% 2*(513-(find(BgAwgn==max(BgAwgn)))) + feAwgn

FdProposedRay = FcoarseRay + FfineRay - 1
FdSchmidlRay = (1+256-find(BgAwgn==max(BgAwgn)))*2 +feAwgn%2*(513-(find(BgRay==max(BgRay)))) + feRay

% figure;plot(abs(conj(F1).*F2));

% figure;plot(abs(FkAwgn));

%%
% figure;hold on;
% plot((fft(A/100)));
% plot((frs),'r');
% figure;plot(abs(RespOfFind));
%%
% figure;plot(Bg);


