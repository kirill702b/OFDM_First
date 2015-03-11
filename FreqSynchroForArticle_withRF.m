%Идеи: 1) не работает для начала потому, что frs2 не используется для
%составления последовательности B, но используется в поиске, заменить
%а потом возможно, это нельзя просто потому что А особенная
%последовательность у которой - что-то типа только четные.
%поэтому можно взять половину B = A(end:-1:1), вторая половина случайная
clear all;close all;
tic;
FreqOffset = 12.4; % Frequency offcet in subcarrier spacing
SNR =30;
%
ProposedError =0;
ShmidlError = 0;
FreqDop = 100;
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
Nnoise = 400;%number of noise samples
phiRand = 2*pi/3;
Nshow = 100;

FullTimingSyncSim = 0;%Enable simulation of second method of timing synchronization
EnableGraphs = 1;
EnableOutput = 1;
% 

frs = (randi(2,[1 N/2])-1.5)*2;%mseq(2,log2(N/2)); последовательность +-1
frs2 = (randi(2,[1 N])-1.5)*2;% вторая последовательность

A = ifft(frs,N/2);%A -последовательность длины N/2
%ifft(frs2,N);% последовательность длины N
AA = [A(1:N/2),A(1:N/2)];% последовательность длины N
B = conj(AA(end:-1:1));
AA_B = [AA, zeros(1,fix(CP*N)),B]; %две последовательности длины 2048+102=2150
seq1 = [zeros(1,Nnoise),real(AA_B),zeros(1,Nnoise)]+1i* [zeros(1,Nnoise),imag(AA_B),zeros(1,Nnoise)];% последовательность с началом и концом без сигнала
%seq1 - это для метода шмидля
%seq2 - это для метода proposed

frs3 = (randi(2,[1 N/4])-1.5)*2;
C = ifft(frs3,N/4);
CCmCmC = [C,C,-C,-C];
CCmCmC_B = [CCmCmC, zeros(1,fix(CP*N)),B];
seq2 = [zeros(1,Nnoise),real(CCmCmC),zeros(1,Nnoise)]+1i* [zeros(1,Nnoise),imag(CCmCmC),zeros(1,Nnoise)];
CiCcCciC=[C,C(end:-1:1),conj(C),conj(C(end:-1:1))];
CiCcCciC_B = [CiCcCciC,zeros(1,fix(CP*N)),B];
seq3 = [zeros(1,Nnoise),real(CiCcCciC),zeros(1,Nnoise)]+1i* [zeros(1,Nnoise),imag(CiCcCciC),zeros(1,Nnoise)];




% длина 6*N+2150 = 8294

seq1Up = resample(seq1,Fs/Fs1,1); %последовательность на частоте 2ГГц, была на 20MHz
seq3Up = resample(seq3,Fs/Fs1,1); 
if FullTimingSyncSim
    seq2Up = resample(seq2,Fs/Fs1,1); 
end
tUp = (1:length(seq1Up))/Fs;% вектор времен для этой новой последовательности 
tUp3 = (1:length(seq3Up))/Fs;
seq1UpFc = seq1Up.*exp(1i*2*pi*(Fc+FreqOffset*Fs1/N)*tUp);% последовательность уже на несущей + доплеровкий сдвиг
seq3UpFc = seq3Up.*exp(1i*2*pi*(Fc+FreqOffset*Fs1/N)*tUp3);
if FullTimingSyncSim
    seq2UpFc = seq2Up.*exp(1i*2*pi*(Fc+FreqOffset*Fs1/N)*tUp3);
end


tt1 = (1:length(AA)*Fs/Fs1)/Fs;
nfft = 2^15;
% [pw,f]=pwelch(seq1UpFcAwgn(3*N*Fs/Fs1+1:3*N*Fs/Fs1+length(AA)*Fs/Fs1),[],[],nfft,Fs);
% figure;plot(f,10*log10(pw));grid on;% 
% [pw,f]=pwelch(seq1UpFcAwgn,[],[],nfft,Fs);
% figure;plot(f,10*log10(pw));grid on;
% [pw,f]=pwelch(aa1,[],[],nfft,Fs);
% figure;plot(f-Fs/2,[pw(nfft/2+1:nfft);pw(1:nfft/2)]);grid on;
if ~exist('RayCh1')
    RayCh1 = rayleighchan(1/Fs,FreqDop,DelOfPath,AvRecPwr);
%     RayCh1 = ricianchan(1/Fs,277.3,10,DelOfPath,AvRecPwr);
    RayCh1.ResetBeforeFiltering = 1;
end
% RayCh1.
seq1UpFcAwgn = awgn(seq1UpFc,SNR);%,'measured');%сигнал после АГБШ канала
seq1UpFcRay = awgn(filter(RayCh1,seq1UpFc),SNR);%,'measured');

seq3UpFcAwgn = awgn(seq3UpFc,SNR,'measured');%сигнал после АГБШ канала
seq3UpFcRay = awgn(filter(RayCh1,seq3UpFc),SNR,'measured');

if FullTimingSyncSim
    seq2UpFcAwgn = awgn(seq2UpFc,SNR,'measured');%сигнал после АГБШ канала
    seq2UpFcRay = awgn(filter(RayCh1,seq2UpFc),SNR,'measured');
end

seq1UpFcAwgnF0 = seq1UpFcAwgn.*exp(-1i*2*pi*Fc*tUp);%сигнал после переноса частоты на 0
seq1UpFcRayF0 = seq1UpFcRay.*exp(-1i*2*pi*Fc*tUp);%сигнал после переноса частоты на 0
seq3UpFcAwgnF0 = seq3UpFcAwgn.*exp(-1i*2*pi*Fc*tUp3);%сигнал после переноса частоты на 0
seq3UpFcRayF0 = seq3UpFcRay.*exp(-1i*2*pi*Fc*tUp3);%сигнал после переноса частоты на 0
if FullTimingSyncSim
    seq2UpFcAwgnF0 = seq2UpFcAwgn.*exp(-1i*2*pi*Fc*tUp3);%сигнал после переноса частоты на 0
    seq2UpFcRayF0 = seq2UpFcRay.*exp(-1i*2*pi*Fc*tUp3);%сигнал после переноса частоты на 0
end
% figure;hold on;grid on;
% plot(seq1UpFcAwgnF0);
% plot(seq1UpFcRayF0,'r');
% [pw,f]=pwelch(seq1UpFcAwgnF0,[],[],nfft,Fs);% спектр сигнала 
% figure;hold on;grid on;plot(f,10*log10(pw));grid on;%рисуем% 
% [pw,f]=pwelch(seq1UpFcRayF0,[],[],nfft,Fs);% спектр сигнала 
% plot(f,10*log10(pw),'r');grid on;%рисуем

seq1UpFcAwgnF0Dec500 = decimate(decimate(decimate(decimate(seq1UpFcAwgnF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
seq1UpFcRayF0Dec500 = decimate(decimate(decimate(decimate(seq1UpFcRayF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
seq1UpFcAwgnF0Dec500Lpf=seq1UpFcAwgnF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcAwgnF0Dec500);%фильтрация для итоговой децимации
seq1UpFcRayF0Dec500Lpf = seq1UpFcRayF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcRayF0Dec500);% фильтрация для итоговой децимации
seq1UpFcAwgnF0Dec500LpfFs1 = decimate(seq1UpFcAwgnF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты
seq1UpFcRayF0Dec500LpfFs1 = decimate(seq1UpFcRayF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты

seq3UpFcAwgnF0Dec500 = decimate(decimate(decimate(decimate(seq3UpFcAwgnF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
seq3UpFcRayF0Dec500 = decimate(decimate(decimate(decimate(seq3UpFcRayF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
seq3UpFcAwgnF0Dec500Lpf=seq3UpFcAwgnF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcAwgnF0Dec500);%фильтрация для итоговой децимации
seq3UpFcRayF0Dec500Lpf = seq3UpFcRayF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcRayF0Dec500);% фильтрация для итоговой децимации
seq3UpFcAwgnF0Dec500LpfFs1 = decimate(seq3UpFcAwgnF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты
seq3UpFcRayF0Dec500LpfFs1 = decimate(seq3UpFcRayF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты



for j=N+1+fix(CP*N)/2:N+2*Nnoise
    P1 = 0;R1 = 0;
    for k = 0:N-1
        P1 = P1 + seq1UpFcAwgnF0Dec500LpfFs1(1,j+k-1+fix(CP*N)/2)*seq1UpFcAwgnF0Dec500LpfFs1(1,j-k-fix(CP*N)/2);
        R1 = R1 + power(abs(seq1UpFcAwgnF0Dec500LpfFs1(1,j+k-1+fix(CP*N)/2)),2);
%         abs(seq3UpFcAwgnF0Dec500LpfFs1(1,j+k-1)*seq3UpFcAwgnF0Dec500LpfFs1(1,j-k));%power(abs(seq3(1,j+k-1)),2);
    end  
    RespOfFind(1,j-N/2-N/2-fix(CP*N)/2) = power(abs(P1),2)/power(R1,2);
end
for j=N+1+fix(CP*N)/2:N+2*Nnoise
    P1 = 0;R1 = 0;
    for k = 0:N-1
        P1 = P1 + seq1UpFcRayF0Dec500LpfFs1(1,j+k-1+fix(CP*N)/2)*seq1UpFcRayF0Dec500LpfFs1(1,j-k-fix(CP*N)/2);
        R1 = R1 + power(abs(seq1UpFcRayF0Dec500LpfFs1(1,j+k-1+fix(CP*N)/2)),2);
%         abs(seq3UpFcAwgnF0Dec500LpfFs1(1,j+k-1)*seq3UpFcAwgnF0Dec500LpfFs1(1,j-k));%power(abs(seq3(1,j+k-1)),2);
    end  
    RespOfFind(2,j-N/2-N/2-fix(CP*N)/2) = power(abs(P1),2)/power(R1,2);
end

for j=N/2+1:N/2+2*Nnoise
    P1 = 0;R1 = 0;
    for k = 0:N/2
        P1 = P1 + seq3UpFcAwgnF0Dec500LpfFs1(1,j+k-1)*seq3UpFcAwgnF0Dec500LpfFs1(1,j-k);
        R1 = R1 + power(abs(seq3UpFcAwgnF0Dec500LpfFs1(1,j+k-1)),2);
%         abs(seq3UpFcAwgnF0Dec500LpfFs1(1,j+k-1)*seq3UpFcAwgnF0Dec500LpfFs1(1,j-k));%power(abs(seq3(1,j+k-1)),2);
    end  
    RespOfFind(3,j-N/2) = power(abs(P1),2)/power(R1,2);
end
for j=N/2+1:N/2+2*Nnoise
    P1 = 0;R1 = 0;
    for k = 0:N/2
        P1 = P1 + seq3UpFcRayF0Dec500LpfFs1(1,j+k-1)*seq3UpFcRayF0Dec500LpfFs1(1,j-k);
        R1 = R1 + power(abs(seq3UpFcRayF0Dec500LpfFs1(1,j+k-1)),2);
%         abs(seq3UpFcRayF0Dec500LpfFs1(1,j+k-1)*seq3UpFcRayF0Dec500LpfFs1(1,j-k));%power(abs(seq3(1,j+k-1)),2);
    end  
    RespOfFind(4,j-N/2) = power(abs(P1),2)/power(R1,2);
end



if FullTimingSyncSim
    seq2UpFcAwgnF0Dec500 = decimate(decimate(decimate(decimate(seq2UpFcAwgnF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
    seq2UpFcRayF0Dec500 = decimate(decimate(decimate(decimate(seq2UpFcRayF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
    seq2UpFcAwgnF0Dec500Lpf=seq2UpFcAwgnF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcAwgnF0Dec500);%фильтрация для итоговой децимации
    seq2UpFcRayF0Dec500Lpf = seq2UpFcRayF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcRayF0Dec500);% фильтрация для итоговой децимации
    seq2UpFcAwgnF0Dec500LpfFs1 = decimate(seq2UpFcAwgnF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты
    seq2UpFcRayF0Dec500LpfFs1 = decimate(seq2UpFcRayF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты
    for j=1:2*Nnoise
        P1 = sum(conj(seq2UpFcAwgnF0Dec500LpfFs1(1,j:j+N/4-1)).*seq2UpFcAwgnF0Dec500LpfFs1(1,j+N/4:j+N/4-1+N/4)) ...
            + sum(conj(seq2UpFcAwgnF0Dec500LpfFs1(1,j+2*N/4:j+N/4-1+2*N/4)).*seq2UpFcAwgnF0Dec500LpfFs1(1,j+N/4+2*N/4:j+N/4-1+N/4+2*N/4));
        R1 = sum(power(abs(seq2UpFcAwgnF0Dec500LpfFs1(1,j+N/4:j+N/4-1+N/4)),2)) + sum(power(abs(seq2UpFcAwgnF0Dec500LpfFs1(1,j+N/4+2*N/4:j+N/4-1+N/4+2*N/4)),2));
%         sum(abs(conj(seq2UpFcAwgnF0Dec500LpfFs1(1,j:j+N/4-1)).*seq2UpFcAwgnF0Dec500LpfFs1(1,j+N/4:j+N/4-1+N/4))) ...
%             + sum(abs(conj(seq2UpFcAwgnF0Dec500LpfFs1(1,j+2*N/4:j+N/4-1+2*N/4)).*seq2UpFcAwgnF0Dec500LpfFs1(1,j+N/4+2*N/4:j+N/4-1+N/4+2*N/4))) ;
        % 
        RespOfFind(5,j) = power(abs(P1),2)/power(R1,2);
    end
    for j=1:2*Nnoise
        P1 = sum(conj(seq2UpFcRayF0Dec500LpfFs1(1,j:j+N/4-1)).*seq2UpFcRayF0Dec500LpfFs1(1,j+N/4:j+N/4-1+N/4)) ...
            + sum(conj(seq2UpFcRayF0Dec500LpfFs1(1,j+2*N/4:j+N/4-1+2*N/4)).*seq2UpFcRayF0Dec500LpfFs1(1,j+N/4+2*N/4:j+N/4-1+N/4+2*N/4));
        R1 = sum(power(abs(seq2UpFcRayF0Dec500LpfFs1(1,j+N/4:j+N/4-1+N/4)),2)) + sum(power(abs(seq2UpFcRayF0Dec500LpfFs1(1,j+N/4+2*N/4:j+N/4-1+N/4+2*N/4)),2));
%         sum(abs(conj(seq2UpFcRayF0Dec500LpfFs1(1,j:j+N/4-1)).*seq2UpFcRayF0Dec500LpfFs1(1,j+N/4:j+N/4-1+N/4))) ...
%             + sum(abs(conj(seq2UpFcRayF0Dec500LpfFs1(1,j+2*N/4:j+N/4-1+2*N/4)).*seq2UpFcRayF0Dec500LpfFs1(1,j+N/4+2*N/4:j+N/4-1+N/4+2*N/4))) ;% 
        RespOfFind(6,j) = power(abs(P1),2)/power(R1,2);
    end
end
tResp = 1:length(RespOfFind);
 close all;
if EnableGraphs
    figure;hold on;grid on;%рисуем отклики во времени
    plot(tResp(:)-Nnoise-1,abs(RespOfFind(1,:)),'--');%2.2*N:end-3*N
    plot(tResp(:)-Nnoise-1,abs(RespOfFind(3,:)),'r');
    if FullTimingSyncSim
        plot(tResp(:)-Nnoise-1,abs(RespOfFind(5,:)),'k:');
        legend('Метод Шмидля','Метод Минна','Метод Парка');
    else
        legend('Метод Шмидля','Метод Минна');
    end
    xlabel('Дискретные отсчеты');ylabel('Временная метрика, ед.');


    figure;hold on;grid on;
    plot(tResp(:)-Nnoise-1,abs(RespOfFind(2,:)),'--');%2.2*N:end-3*N
    plot(tResp(:)-Nnoise-1,abs(RespOfFind(4,:)),'r');
    if FullTimingSyncSim
        plot(tResp(:)-Nnoise-1,abs(RespOfFind(6,:)),'k:');
        legend('Метод Шмидля','Метод Минна','Метод Парка');
    else
        legend('Метод Шмидля','Метод Минна');
    end
    xlabel('Дискретные отсчеты');ylabel('Временная метрика, ед.');
    drawnow;
end
maxAwgnSch = 0;
maxAwgnIndSch = 0;
maxRaySch = 0;
maxRayIndSch = 0;
curMax = 0;
for j=1:length(tResp)-103
    curMax = sum(RespOfFind(1,j:j+102));
    if curMax>maxAwgnSch
        maxAwgnSch = curMax;
        maxAwgnIndSch = j;
    end
    curMax = sum(RespOfFind(2,j:j+102));
    if curMax>maxRaySch
        maxRaySch = curMax;
        maxRayIndSch = j;
    end
end
maxAwgnIndPro = find(RespOfFind(3,:)==max(RespOfFind(3,:)));
maxRayIndPro = find(RespOfFind(4,:)==max(RespOfFind(4,:)));

if EnableOutput
    maxAwgnIndSch-Nnoise-1
    maxRayIndSch-Nnoise-1
    maxAwgnIndPro-Nnoise-1
    maxRayIndPro-Nnoise-1
end
if maxAwgnIndPro-Nnoise-1~=0 || maxRayIndPro-Nnoise-1~=0
    maxAwgnIndPro
    maxRayIndPro
    ProposedError = 1;
end
% 

Vk = sqrt(2).*frs2(1:2:end)./frs;%вспомогательная последовательность
% figure;plot(abs(seq1));
t1 = (1:length(AA(1,:)));
t12 = (1:length(AA_B(1,:)));

% seq1df(1,:) = seq1(1,3*N+1:3*N+N).*exp(1i*(2*pi*FreqOffset*t1/N));%Сначала брал только неискаженную последовательность, сдвинутую
% seq12df(1,:) = seq1(1,3*N+1:3*N+N+N+fix(CP*N)).*exp(1i*(2*pi*FreqOffset*t12/N));

%%
seq1df(1,:) = seq1UpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndPro:maxAwgnIndPro+N-1);%Теперь уже все по честному 
seq12df(1,:) = seq1UpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndPro:maxAwgnIndPro-1+N+N+fix(CP*N));%АГБШ

seq1df(2,:) = seq1UpFcRayF0Dec500LpfFs1(1,maxRayIndPro:maxRayIndPro-1+N);%Теперь уже все по честному
seq12df(2,:) = seq1UpFcRayF0Dec500LpfFs1(1,maxRayIndPro:maxRayIndPro-1+N+N+fix(CP*N));%РЭЛЕЙ

% FkAwgn = fft(conj(AA).*seq1df(1,:),N);
% FkAwgn = fft(conj(B).*seq12df(1,N+fix(CP*N)+1:N+fix(CP*N)+N),N);
FkAwgn = fft(conj(AA).*seq12df(1,1:N),N);
FkAwgn = [FkAwgn(N/2+1:end),FkAwgn(1:N/2)];
kmax =  find(abs(FkAwgn)==max(abs(FkAwgn)));
FcoarseAwgn = kmax-N/2;
if abs(FkAwgn(kmax-1))<= abs(FkAwgn(kmax+1))
    alp = 1;
else
    alp = -1;
end
FfineAwgn = alp/(abs(FkAwgn(kmax))/abs(FkAwgn(kmax+alp))+1);


% FkRay = fft(conj(AA).*seq1df(2,:),N);
% FkRay = fft(conj(B).*seq12df(2,N+fix(CP*N)+1:N+fix(CP*N)+N),N);
FkRay = fft(conj(AA).*seq12df(2,1:N),N);
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

seq1df(1,:) = seq1UpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndSch:maxAwgnIndSch+N-1);%Теперь уже все по честному 
seq12df(1,:) = seq1UpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndSch:maxAwgnIndSch-1+N+N+fix(CP*N));%АГБШ

seq1df(2,:) = seq1UpFcRayF0Dec500LpfFs1(1,maxRayIndSch:maxRayIndSch-1+N);%Теперь уже все по честному
seq12df(2,:) = seq1UpFcRayF0Dec500LpfFs1(1,maxRayIndSch:maxRayIndSch-1+N+N+fix(CP*N));%РЭЛЕЙ
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
FdProposedAwgn = FcoarseAwgn + FfineAwgn - 1;
FdSchmidlAwgn = (1+256-find(BgAwgn==max(BgAwgn)))*2 +feAwgn;% 2*(513-(find(BgAwgn==max(BgAwgn)))) + feAwgn

FdProposedRay = FcoarseRay + FfineRay - 1;
FdSchmidlRay = (1+256-find(BgRay==max(BgRay)))*2 +feRay;%2*(513-(find(BgRay==max(BgRay)))) + feRay
if abs(FdSchmidlRay-FreqOffset )>0.9
    ShmidlError = 1;
end
if EnableOutput
   fprintf('Error of proposed method is %.5f (awgn), %.5f (rayleygh)\n',FdProposedAwgn-FreqOffset,FdProposedRay-FreqOffset);
   fprintf('Error of Schmidl method is %.5f (awgn), %.5f (rayleygh)\n',FdSchmidlAwgn-FreqOffset,FdSchmidlRay-FreqOffset);      
   
end
% toc;
% figure;plot(abs(conj(F1).*F2));

% figure;plot(abs(FkAwgn));

%%
% figure;hold on;
% plot((fft(A/100)));
% plot((frs),'r');
% figure;plot(abs(RespOfFind));
%%
% figure;plot(Bg);


