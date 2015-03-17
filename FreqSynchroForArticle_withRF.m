%Почему бы не провести исследование зависимости от частоты доплера, от
%часотного сдвига, параметров канала рэлея. Причем не нужно хвататься за
%выбросы. просто усредняю и все
clear all;close all;
FreqOffset = 12.4; % Frequency offcet in subcarrier spacing
SNR = 30;
FreqDop = 27;
%
% tic;
Nfft=2^17;
dev = 50;
EnableGraphs = 0;
EnableOutput = 0;
PartOfB1 = 0.0;
ProposedError =0;
SchmidlError = 0;
SchmidlOldError = 0;
MaxFreqError = 0.5;

Nc = 1000; % Number of subcarriers
N = 1024; % Number of IFFT points
Rd = 18e6; % Data rate bit per second
CP = 0.1; % Cyclic prefix length as part of symbol period
Fs1 =(Rd/(1-CP));%20mHz, 50ns
% 
NumOfOsc = 16; % Number of Oscillator for Jakes' model
Fc = 2e8; % Hz, Carrier frequency 2GHz
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

%
frs3 = (randi(2,[1 N])-1.5)*2;% вторая последовательность
frs4 = (randi(2,[1 N/4])-1.5)*2;
C = ifft(frs4,N/4);
BP = ifft(frs3,N);

CiCcCciC=[C,C(end:-1:1),conj(C),conj(C(end:-1:1))];
CiCcCciC_B = [CiCcCciC,zeros(1,fix(CP*N)),BP];
frs = (randi(2,[1 N/2])-1.5)*2;%mseq(2,log2(N/2)); последовательность +-1

frs2Old = (randi(2,[1 N])-1.5)*2;
G1 = ifft(frs3,N);
% G2 = ifft((randi(2,[1 N/2])-1.5)*2);
G2 =G1(1:N/2);
% GciG=[G2,G2(end:-1:1),zeros(1,fix(CP*N)),conj(G2),conj(G2(end:-1:1))];
GciG=CiCcCciC_B;%[CiCcCciC,zeros(1,fix(CP*N)),CiCcCciC];
G3=GciG(end-N+1:end);
seq1 = [zeros(1,Nnoise),real(GciG),zeros(1,Nnoise)]+1i* [zeros(1,Nnoise),imag(GciG),zeros(1,Nnoise)];
A = ifft(frs,N/2);%A -последовательность длины N/2
%frs3(1:(N-Nc)/2) = zeros(1,(N-Nc)/2);frs3(end:-1:end-(N-Nc)/2+1) = zeros(1,(N-Nc)/2);
B1 = ifft(frs3,N);% последовательность длины N
AA = [A(1:N/2),A(1:N/2)];% последовательность длины N
% B = [conj(A(end:-1:1)), B1(1:N/2)];
% B = [conj(AA(end:-1:1+fix(N*PartOfB1))), B1(1:fix(N*PartOfB1))/PartOfB1];%
B = ifft(frs3,N);
BOld = ifft(frs2Old,N);
frs2 = fft(B,N);
AA_B = [AA, zeros(1,fix(CP*N)),B]; %две последовательности длины 2048+102=2150
% seq1 = [zeros(1,Nnoise),real(AA_B),zeros(1,Nnoise)]+1i* [zeros(1,Nnoise),imag(AA_B),zeros(1,Nnoise)];% последовательность с началом и концом без сигнала
AA_BOld = [AA, zeros(1,fix(CP*N)),BOld];
seq1Old = [zeros(1,Nnoise),real(AA_BOld),zeros(1,Nnoise)]+1i* [zeros(1,Nnoise),imag(AA_BOld),zeros(1,Nnoise)];% последовательность с началом и концом без сигнала
%seq1 - это для метода шмидля
%seq2 - это для метода proposed


CCmCmC = [C,C,-C,-C];
CCmCmC_B = [CCmCmC, zeros(1,fix(CP*N)),BP];
seq2 = [zeros(1,Nnoise),real(CCmCmC),zeros(1,Nnoise)]+1i* [zeros(1,Nnoise),imag(CCmCmC),zeros(1,Nnoise)];

seq3 = [zeros(1,Nnoise),real(CiCcCciC_B),zeros(1,Nnoise)]+1i* [zeros(1,Nnoise),imag(CiCcCciC_B),zeros(1,Nnoise)];




% длина 6*N+2150 = 8294

seq1Up = resample(seq1,Fs/Fs1,1); %последовательность на частоте 2ГГц, была на 20MHz
seq1OldUp = resample(seq1Old,Fs/Fs1,1);
seq3Up = resample(seq3,Fs/Fs1,1); 
if FullTimingSyncSim
    seq2Up = resample(seq2,Fs/Fs1,1); 
end
tUp = (1:length(seq1Up))/Fs;% вектор времен для этой новой последовательности 
tUp3 = (1:length(seq3Up))/Fs;
seq1UpFc = seq1Up.*exp(1i*2*pi*(Fc+FreqOffset*Fs1/N)*tUp);% последовательность уже на несущей + доплеровкий сдвиг
seq1OldUpFc = seq1OldUp.*exp(1i*2*pi*(Fc+FreqOffset*Fs1/N)*tUp);% последовательность уже на несущей + доплеровкий сдвиг
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
% RayCh1
seq1UpFcAwgn = awgn(seq1UpFc,SNR);%,'measured');%сигнал после АГБШ канала
seq1UpFcRay = awgn(filter(RayCh1,seq1UpFc),SNR);%,'measured');

seq1OldUpFcAwgn = awgn(seq1OldUpFc,SNR);%,'measured');%сигнал после АГБШ канала
seq1OldUpFcRay = awgn(filter(RayCh1,seq1OldUpFc),SNR);%,'measured');

seq3UpFcAwgn = awgn(seq3UpFc,SNR);%сигнал после АГБШ канала
seq3UpFcRay = awgn(filter(RayCh1,seq3UpFc),SNR);

if FullTimingSyncSim
    seq2UpFcAwgn = awgn(seq2UpFc,SNR);%сигнал после АГБШ канала
    seq2UpFcRay = awgn(filter(RayCh1,seq2UpFc),SNR);
end

seq1UpFcAwgnF0 = seq1UpFcAwgn.*exp(-1i*2*pi*Fc*tUp);%сигнал после переноса частоты на 0
seq1UpFcRayF0 = seq1UpFcRay.*exp(-1i*2*pi*Fc*tUp);%сигнал после переноса частоты на 0
seq1OldUpFcAwgnF0 = seq1OldUpFcAwgn.*exp(-1i*2*pi*Fc*tUp);%сигнал после переноса частоты на 0
seq1OldUpFcRayF0 = seq1OldUpFcRay.*exp(-1i*2*pi*Fc*tUp);%сигнал после переноса частоты на 0

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

Dec = Fs/Fs1;
if Dec-fix(Dec)~=0
    display('Decimation factor is not integer!!! Exit.');
    return;
end
Decs = sort(factor(Dec),'descend');


seq1UpFcAwgnF0Dec500 = decimate(seq1UpFcAwgnF0,Decs(1),64,'fir');% децимация в 500раз (всего было 1000)
seq1UpFcRayF0Dec500 = decimate(seq1UpFcRayF0,Decs(1),64,'fir');% децимация в 500раз (всего было 1000)
seq1OldUpFcAwgnF0Dec500 = decimate(seq1OldUpFcAwgnF0,Decs(1),64,'fir');% децимация в 500раз (всего было 1000)
seq1OldUpFcRayF0Dec500 = decimate(seq1OldUpFcRayF0,Decs(1),64,'fir');% децимация в 500раз (всего было 1000)
seq3UpFcAwgnF0Dec500 = decimate(seq3UpFcAwgnF0,Decs(1),64,'fir');% децимация в 500раз (всего было 1000)
seq3UpFcRayF0Dec500 = decimate(seq3UpFcRayF0,Decs(1),64,'fir');% децимация в 500раз (всего было 1000)
for j=2:length(Decs)-1
    seq1UpFcAwgnF0Dec500 = decimate(seq1UpFcAwgnF0Dec500,Decs(j),64,'fir');% децимация в 500раз (всего было 1000)
    seq1UpFcRayF0Dec500 = decimate(seq1UpFcRayF0Dec500,Decs(j),64,'fir');% децимация в 500раз (всего было 1000)
    seq1OldUpFcAwgnF0Dec500 = decimate(seq1OldUpFcAwgnF0Dec500,Decs(j),64,'fir');% децимация в 500раз (всего было 1000)
    seq1OldUpFcRayF0Dec500 = decimate(seq1OldUpFcRayF0Dec500,Decs(j),64,'fir');% децимация в 500раз (всего было 1000)
    seq3UpFcAwgnF0Dec500 = decimate(seq3UpFcAwgnF0Dec500,Decs(j),64,'fir');% децимация в 500раз (всего было 1000)
    seq3UpFcRayF0Dec500 = decimate(seq3UpFcRayF0Dec500,Decs(j),64,'fir');% децимация в 500раз (всего было 1000)
end
% seq1UpFcAwgnF0Dec500 = decimate(decimate(decimate(decimate(seq1UpFcAwgnF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
% seq1UpFcRayF0Dec500 = decimate(decimate(decimate(decimate(seq1UpFcRayF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
% seq1OldUpFcAwgnF0Dec500 = decimate(decimate(decimate(decimate(seq1OldUpFcAwgnF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
% seq1OldUpFcRayF0Dec500 = decimate(decimate(decimate(decimate(seq1OldUpFcRayF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)

if (Decs(end)~=2)
    display('Last decimation factor is not 2!!! Exit.');
    return;
end
seq1UpFcAwgnF0Dec500Lpf=seq1UpFcAwgnF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcAwgnF0Dec500);%фильтрация для итоговой децимации
seq1UpFcRayF0Dec500Lpf = seq1UpFcRayF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcRayF0Dec500);% фильтрация для итоговой децимации
seq1OldUpFcAwgnF0Dec500Lpf=seq1OldUpFcAwgnF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcAwgnF0Dec500);%фильтрация для итоговой децимации
seq1OldUpFcRayF0Dec500Lpf = seq1OldUpFcRayF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcRayF0Dec500);% фильтрация для итоговой децимации

seq1UpFcAwgnF0Dec500LpfFs1 = decimate(seq1UpFcAwgnF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты
seq1UpFcRayF0Dec500LpfFs1 = decimate(seq1UpFcRayF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты
seq1OldUpFcAwgnF0Dec500LpfFs1 = decimate(seq1OldUpFcAwgnF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты
seq1OldUpFcRayF0Dec500LpfFs1 = decimate(seq1OldUpFcRayF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты

% seq3UpFcAwgnF0Dec500 = decimate(decimate(decimate(decimate(seq3UpFcAwgnF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
% seq3UpFcRayF0Dec500 = decimate(decimate(decimate(decimate(seq3UpFcRayF0,5,64,'fir'),5,64,'fir'),5,64,'fir'),4,64,'fir');% децимация в 500раз (всего было 1000)
seq3UpFcAwgnF0Dec500Lpf=seq3UpFcAwgnF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcAwgnF0Dec500);%фильтрация для итоговой децимации
seq3UpFcRayF0Dec500Lpf = seq3UpFcRayF0Dec500;%LPF_fs40MHz_20Mpass_21Mstop(seq1UpFcRayF0Dec500);% фильтрация для итоговой децимации
seq3UpFcAwgnF0Dec500LpfFs1 = decimate(seq3UpFcAwgnF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты
seq3UpFcRayF0Dec500LpfFs1 = decimate(seq3UpFcRayF0Dec500Lpf,2,128,'fir'); % итоговая децимация в 2раза, до изначальной частоты


% RespOfFind:
% 1 - Schmild new AWGN, 2 - Schmidl new Rayleygh
% 3 - Park AWGN, 4 - Park Rayleygh
% 5 - Minn AWGN, 6 - Minn Rayleygh
% 7 - Schmild old AWGN, 8 - Schmidl old Rayleygh


% for j=N/2+1:N/2+2*Nnoise
%     P1 = 0;R1 = 0;
%     for k = 0:N/2
%         P1 = P1 + seq1UpFcAwgnF0Dec500LpfFs1(1,j+k-1)*seq1UpFcAwgnF0Dec500LpfFs1(1,j-k)+ seq1UpFcAwgnF0Dec500LpfFs1(1,j+k-1+N+fix(CP*N))*seq1UpFcAwgnF0Dec500LpfFs1(1,j-k+N+fix(CP*N));
%         R1 = R1 + power(abs(seq1UpFcAwgnF0Dec500LpfFs1(1,j+k-1)),2)+ power(abs(seq1UpFcAwgnF0Dec500LpfFs1(1,j+k-1+N+fix(CP*N))),2);
% %         abs(seq3UpFcAwgnF0Dec500LpfFs1(1,j+k-1)*seq3UpFcAwgnF0Dec500LpfFs1(1,j-k));%power(abs(seq3(1,j+k-1)),2);
%     end  
%     RespOfFind(1,j-N/2) = power(abs(P1),2)/power(R1,2);
% end
% for j=N/2+1:N/2+2*Nnoise
%     P1 = 0;R1 = 0;
%     for k = 0:N/2
%         P1 = P1 + seq1UpFcRayF0Dec500LpfFs1(1,j+k-1)*seq1UpFcRayF0Dec500LpfFs1(1,j-k)+ seq1UpFcRayF0Dec500LpfFs1(1,j+k-1+N+fix(CP*N))*seq1UpFcRayF0Dec500LpfFs1(1,j-k+N+fix(CP*N));
%         R1 = R1 + power(abs(seq1UpFcRayF0Dec500LpfFs1(1,j+k-1)),2)+ power(abs(seq1UpFcRayF0Dec500LpfFs1(1,j+k-1+N+fix(CP*N))),2);
% %         abs(seq3UpFcRayF0Dec500LpfFs1(1,j+k-1)*seq3UpFcRayF0Dec500LpfFs1(1,j-k));%power(abs(seq3(1,j+k-1)),2);
%     end  
%     RespOfFind(2,j-N/2) = power(abs(P1),2)/power(R1,2);
% end

for j=N/2+1:N/2+2*Nnoise
    P1 = 0;R1 = 0;
    for k = 0:N/2
        P1 = P1 + seq1UpFcAwgnF0Dec500LpfFs1(1,j+k-1)*seq1UpFcAwgnF0Dec500LpfFs1(1,j-k);
        R1 = R1 + power(abs(seq1UpFcAwgnF0Dec500LpfFs1(1,j+k-1)),2);
%         abs(seq3UpFcAwgnF0Dec500LpfFs1(1,j+k-1)*seq3UpFcAwgnF0Dec500LpfFs1(1,j-k));%power(abs(seq3(1,j+k-1)),2);
    end  
    RespOfFind(1,j-N/2) = power(abs(P1),2)/power(R1,2);
end
for j=N/2+1:N/2+2*Nnoise
    P1 = 0;R1 = 0;
    for k = 0:N/2
        P1 = P1 + seq1UpFcRayF0Dec500LpfFs1(1,j+k-1)*seq1UpFcRayF0Dec500LpfFs1(1,j-k);
        R1 = R1 + power(abs(seq1UpFcRayF0Dec500LpfFs1(1,j+k-1)),2);
%         abs(seq3UpFcRayF0Dec500LpfFs1(1,j+k-1)*seq3UpFcRayF0Dec500LpfFs1(1,j-k));%power(abs(seq3(1,j+k-1)),2);
    end  
    RespOfFind(2,j-N/2) = power(abs(P1),2)/power(R1,2);
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
%THE SAME SEQUENCE FOR SYNCHRO BOTH SCHMIDL!!!
for j=1:2*Nnoise %???? ?????? (?????????????) ??? ???? ??????
    P1 = sum(conj(seq1OldUpFcAwgnF0Dec500LpfFs1(1,j:j+N/2-1)).*seq1OldUpFcAwgnF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2));
    R1 = sum(power(abs(seq1OldUpFcAwgnF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2)),2));
    RespOfFind(7,j) = power(abs(P1),2)/power(R1,2); % ?????? ?? ?????
end
for j=1:2*Nnoise %???? ?????? (?????????????) ??? ???? ??????
    P1 = sum(conj(seq1OldUpFcRayF0Dec500LpfFs1(1,j:j+N/2-1)).*seq1OldUpFcRayF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2));
    R1 = sum(power(abs(seq1OldUpFcRayF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2)),2));
    RespOfFind(8,j) = power(abs(P1),2)/power(R1,2); % ?????? ?? ?????
end
% for j=1:2*Nnoise %???? ?????? (?????????????) ??? ???? ??????
%     P1 = sum(conj(seq1OldUpFcAwgnF0Dec500LpfFs1(1,j:j+N/2-1)).*seq1OldUpFcAwgnF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2));
%     R1 = sum(power(abs(seq1OldUpFcAwgnF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2)),2));
%     RespOfFind(7,j) = power(abs(P1),2)/power(R1,2); % ?????? ?? ?????
% end
% for j=1:2*Nnoise %???? ?????? (?????????????) ??? ???? ??????
%     P1 = sum(conj(seq1OldUpFcRayF0Dec500LpfFs1(1,j:j+N/2-1)).*seq1OldUpFcRayF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2));
%     R1 = sum(power(abs(seq1OldUpFcRayF0Dec500LpfFs1(1,j+N/2:j+N/2-1+N/2)),2));
%     RespOfFind(8,j) = power(abs(P1),2)/power(R1,2); % ?????? ?? ?????
% end


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
    plot(tResp(:)-Nnoise-1,abs(RespOfFind(7,:)),'g');
    if FullTimingSyncSim
        plot(tResp(:)-Nnoise-1,abs(RespOfFind(5,:)),'k:');
        legend('Метод Минна для Шмидля','Метод Минна','Метод Шмидля','Метод Парка');
    else
        legend('Метод Минна для Шмидля','Метод Минна','Метод Шмидля');
    end
    xlabel('Дискретные отсчеты');ylabel('Временная метрика в канале с АГБШ, ед.');


    figure;hold on;grid on;
    plot(tResp(:)-Nnoise-1,abs(RespOfFind(2,:)),'--');%2.2*N:end-3*N
    plot(tResp(:)-Nnoise-1,abs(RespOfFind(4,:)),'r');
    plot(tResp(:)-Nnoise-1,abs(RespOfFind(8,:)),'g');
    if FullTimingSyncSim
        plot(tResp(:)-Nnoise-1,abs(RespOfFind(6,:)),'k:');
        legend('Метод Минна для Шмидля','Метод Минна','Метод Шмидля','Метод Парка');
    else
        legend('Метод Минна для Шмидля','Метод Минна','Метод Шмидля');
    end
    xlabel('Дискретные отсчеты');ylabel('Временная метрика в Рэлеевском канале, ед.');
    drawnow;
end
%
maxAwgnSchOld = 0;
maxAwgnIndSchOld = 0;
maxRaySchOld = 0;
maxRayIndSchOld = 0;
curMaxOld = 0;
for j=1:length(tResp)-103
    curMaxOld = sum(RespOfFind(7,j:j+102));
    if curMaxOld>maxAwgnSchOld
        maxAwgnSchOld = curMaxOld;
        maxAwgnIndSchOld = j;
    end
    curMaxOld = sum(RespOfFind(8,j:j+102));
    if curMaxOld>maxRaySchOld
        maxRaySchOld = curMaxOld;
        maxRayIndSchOld = j;
    end
end
%
maxAwgnSch = 0;
maxAwgnIndSch = 0;
maxRaySch = 0;
maxRayIndSch = 0;
curMax = 0;
maxAwgnIndSch = find(RespOfFind(1,:)==max(RespOfFind(1,:)));
maxRayIndSch = find(RespOfFind(2,:)==max(RespOfFind(2,:)));
maxAwgnIndPro = find(RespOfFind(3,:)==max(RespOfFind(3,:)));
maxRayIndPro = find(RespOfFind(4,:)==max(RespOfFind(4,:)));

if EnableOutput
    fprintf('Time synchronization error Schmidl new %d %d, proposed %d %d, Schmidl old %d %d\n',maxAwgnIndSch-Nnoise-1,maxRayIndSch-Nnoise-1,maxAwgnIndPro-Nnoise-1,...
        maxRayIndPro-Nnoise-1,maxAwgnIndSchOld-Nnoise-1,maxRayIndSchOld-Nnoise-1);   
end
if maxAwgnIndPro-Nnoise-1~=0 || maxRayIndPro-Nnoise-1~=0
    fprintf('ERROR, Proposed!!!\n ');
%     maxAwgnIndPro
%     maxRayIndPro
    ProposedError = 1;
end

% 

% Vk = sqrt(2).*frs2(1:2:end)./frs;%вспомогательная последовательность
% VkOld = sqrt(2).*frs2Old(1:2:end)./frs;
% Vk1 = frs2(1:2:end)./frs;%вспомогательная последовательность
% Vk2 = frs2(2:2:end)./frs;
VkOld = frs2Old(1:2:end)./frs;
% figure;plot(abs(seq1));
t1 = (1:length(AA(1,:)));
t12 = (1:length(AA_B(1,:)));

% seq1df(1,:) = seq1(1,3*N+1:3*N+N).*exp(1i*(2*pi*FreqOffset*t1/N));%Сначала брал только неискаженную последовательность, сдвинутую
% seq12df(1,:) = seq1(1,3*N+1:3*N+N+N+fix(CP*N)).*exp(1i*(2*pi*FreqOffset*t12/N));

%% Proposed
seq1df(1,:) = seq3UpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndPro:maxAwgnIndPro+N-1);%Теперь уже все по честному 
seq12df(1,:) = seq3UpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndPro:maxAwgnIndPro-1+N+N+fix(CP*N));%АГБШ

seq1df(2,:) = seq3UpFcRayF0Dec500LpfFs1(1,maxRayIndPro:maxRayIndPro-1+N);%Теперь уже все по честному
seq12df(2,:) = seq3UpFcRayF0Dec500LpfFs1(1,maxRayIndPro:maxRayIndPro-1+N+N+fix(CP*N));%РЭЛЕЙ

% FkAwgn = fft(conj(AA).*seq1df(1,:),N);
% FkAwgn = fft(conj(B).*seq12df(1,N+fix(CP*N)+1:N+fix(CP*N)+N),N);
FkAwgn = fft(conj(BP).*seq12df(1,1+N+fix(CP*N):N+N+fix(CP*N)),N);
FkAwgn = [FkAwgn(N/2+1:end),FkAwgn(1:N/2)];
kmax =  find(abs(FkAwgn)==max(abs(FkAwgn)));
FcoarseAwgn = kmax-N/2;
if kmax>=2 && kmax <=1023
    if abs(FkAwgn(kmax-1))<= abs(FkAwgn(kmax+1))
        alp = 1;
    else
        alp = -1;
    end
    FfineAwgn = alp/(abs(FkAwgn(kmax))/abs(FkAwgn(kmax+alp))+1);
end

% FkRay = fft(conj(AA).*seq1df(2,:),N);
% FkRay = fft(conj(B).*seq12df(2,N+fix(CP*N)+1:N+fix(CP*N)+N),N);
FkRay = fft(conj(BP).*seq12df(2,1+N+fix(CP*N):N+N+fix(CP*N)),N);
FkRay = [FkRay(N/2+1:end),FkRay(1:N/2)];
kmax =  find(abs(FkRay)==max(abs(FkRay)));
FcoarseRay = kmax-N/2;
if kmax>=2 && kmax <=1023
    if abs(FkRay(kmax-1))<= abs(FkRay(kmax+1))
        alp = 1;
    else
        alp = -1;
    end
    FfineRay = alp/(abs(FkRay(kmax))/abs(FkRay(kmax+alp))+1);
end    
FdProposedAwgn = FcoarseAwgn + FfineAwgn - 1;
FdProposedRay = FcoarseRay + FfineRay - 1;



% Schmidl new ___________________________________________________
seq1df(1,:) = seq1UpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndSch:maxAwgnIndSch+N-1);%Теперь уже все по честному 
seq12df(1,:) = seq1UpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndSch:maxAwgnIndSch-1+N+N+fix(CP*N));%АГБШ

seq1df(2,:) = seq1UpFcRayF0Dec500LpfFs1(1,maxRayIndSch:maxRayIndSch-1+N);%Теперь уже все по честному
seq12df(2,:) = seq1UpFcRayF0Dec500LpfFs1(1,maxRayIndSch:maxRayIndSch-1+N+N+fix(CP*N));%РЭЛЕЙ



FkAwgn = fft(conj(G3).*seq12df(1,1+N+fix(CP*N):N+N+fix(CP*N)),N);
FkAwgn = [FkAwgn(N/2+1:end),FkAwgn(1:N/2)];
kmax =  find(abs(FkAwgn)==max(abs(FkAwgn)));


FcoarseAwgn = kmax-N/2;
if kmax>=2 && kmax <=1023
    if abs(FkAwgn(kmax-1))<= abs(FkAwgn(kmax+1))
        alp = 1;
    else
        alp = -1;
    end
    FfineAwgn = alp/(abs(FkAwgn(kmax))/abs(FkAwgn(kmax+alp))+1);
end

FdSchmidlAwgn = FcoarseAwgn + FfineAwgn - 1;
fmax1 = ((FcoarseAwgn-1-alp)/N*Nfft);
fmax2 =  ((FcoarseAwgn-1+alp)/N*Nfft);
kk = 1;
ts = conj(G3).*seq12df(1,1+N+fix(CP*N):N+N+fix(CP*N)).*window(@hamming,N).';%hamming
for jj = min(fmax1,fmax2)-dev:max(fmax1,fmax2)+dev
    FHRAll(kk) = sum(ts.*exp(-1i*2*pi/Nfft*jj*(0:N-1)));
    kk=kk+1;
end
FHR = FHRAll(dev+1:end-dev);
indm = find(abs(FHR(1:end)) == max(abs(FHR)))+dev;
indmm = sum(sqrt(abs(FHRAll(indm-dev:indm+dev))).*(indm-dev:indm+dev)) / sum(sqrt(abs(FHRAll(indm-dev:indm+dev))));
FdSchmidlAwgn = FcoarseAwgn-2 + (indmm-1 -dev)*N/Nfft;


FkRay = fft(conj(G3).*seq12df(2,1+N+fix(CP*N):N+N+fix(CP*N)),N);
FkRay = [FkRay(N/2+1:end),FkRay(1:N/2)];
kmax =  find(abs(FkRay)==max(abs(FkRay)));
FcoarseRay = kmax-N/2;
if kmax>=2 && kmax <=1023
    if abs(FkRay(kmax-1))<= abs(FkRay(kmax+1))
        alp = 1;
    else
        alp = -1;
    end
    FfineRay = alp/(abs(FkRay(kmax))/abs(FkRay(kmax+alp))+1);
end    

FdSchmidlRay = FcoarseRay + FfineRay - 1;
fmax1 = ((FcoarseRay-1-alp)/N*Nfft);
fmax2 =  ((FcoarseRay-1+alp)/N*Nfft);
kk = 1;
ts = conj(G3).*seq12df(2,1+N+fix(CP*N):N+N+fix(CP*N)).*window(@hamming,N).';
for jj = min(fmax1,fmax2)-dev:max(fmax1,fmax2)+dev
    FHRAll(kk) = sum(ts.*exp(-1i*2*pi/Nfft*jj*(0:N-1)));
    kk=kk+1;
end
FHR = FHRAll(dev+1:end-dev);
indm = find(abs(FHR(1:end)) == max(abs(FHR)))+dev;
indmm = sum(abs(FHRAll(indm-dev:indm+dev)).*(indm-dev:indm+dev)) / sum(abs(FHRAll(indm-dev:indm+dev)));
FdSchmidlRay = FcoarseRay-2 + (indmm-1-dev)*N/Nfft ;
% 
% 
% P1 = sum(conj(seq1df(1,1:N/2)).*seq1df(1,1+N/2:N));
% phi = angle(P1);
% feAwgn = phi/pi;
% Corr = exp(-1i*2*phi*t12/N);
% seq12dfCorr = seq12df(1,:).*Corr;
% F1 = fft(seq12dfCorr(1:N));F1 = [F1(N/2+1:end),F1(1:N/2)];
% F2 = fft(seq12dfCorr(N+fix(CP*N)+1:end));F2 = [F2(N/2+1:end),F2(1:N/2)];
% F1F2 = conj(F1).*[ F2(2:end),F2(1) ];
% F1F2first = conj(F1).*F2;
% % B(1) = power(abs(sum(F1F2.*conj(frs(1:2:end)))),2)/2/power(sum(power(F2,2)),2);
% for i = 1:N/2
%     F1F2sh = circshift(F1F2,[0,2*(i-1)]);
%     F1F2shfirst = circshift(F1F2first,[0,2*(i-1)]);
%     BgRay(i) = power(abs(sum(  F1F2sh(1:2:end).*  conj(Vk2)      )),2)/2/power(sum(power(abs(F2),2)),2)+...
%         power(abs(sum(  F1F2shfirst(1:2:end).*  conj(Vk1)      )),2)/2/power(sum(power(abs(F2),2)),2);
% end
% P1 = sum(conj(seq1df(2,1:N/2)).*seq1df(2,1+N/2:N));
% phi = angle(P1);
% feRay = phi/pi;
% Corr = exp(-1i*2*phi*t12/N);
% seq12dfCorr = seq12df(2,:).*Corr;
% F1 = fft(seq12dfCorr(1:N));F1 = [F1(N/2+1:end),F1(1:N/2)];
% F2 = fft(seq12dfCorr(N+fix(CP*N)+1:end));F2 = [F2(N/2+1:end),F2(1:N/2)];
% F1F2 = conj(F1).*[ F2(2:end),F2(1) ];
% F1F2first = conj(F1).*F2;
% % B(1) = power(abs(sum(F1F2.*conj(frs(1:2:end)))),2)/2/power(sum(power(F2,2)),2);
% for i = 1:N/2
%     F1F2sh = circshift(F1F2,[0,2*(i-1)]);
%     F1F2shfirst = circshift(F1F2first,[0,2*(i-1)]);
%     BgAwgn(i) = power(abs(sum(  F1F2sh(1:2:end).*  conj(Vk2)      )),2)/2/power(sum(power(abs(F2),2)),2)+...
%         power(abs(sum(  F1F2shfirst(1:2:end).*  conj(Vk1)      )),2)/2/power(sum(power(abs(F2),2)),2);
% end
% FdSchmidlAwgn = (1+256-find(BgAwgn==max(BgAwgn)))*2 +feAwgn;% 2*(513-(find(BgAwgn==max(BgAwgn)))) + feAwgn
% 
% FdSchmidlRay = (1+256-find(BgRay==max(BgRay)))*2 +feRay;%2*(513-(find(BgRay==max(BgRay)))) + feRay
if abs(FdSchmidlRay-FreqOffset )>MaxFreqError
    SchmidlError = 1;
    fprintf('ERROR, Schmidl New!!!\n');
end
%Schmidl old _______________________________________________________
seq1Olddf(1,:) = seq1OldUpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndSchOld:maxAwgnIndSchOld+N-1);%Теперь уже все по честному 
seq12Olddf(1,:) = seq1OldUpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndSchOld:maxAwgnIndSchOld-1+N+N+fix(CP*N));%АГБШ

seq1Olddf(2,:) = seq1OldUpFcRayF0Dec500LpfFs1(1,maxRayIndSchOld:maxRayIndSchOld-1+N);%Теперь уже все по честному
seq12Olddf(2,:) = seq1OldUpFcRayF0Dec500LpfFs1(1,maxRayIndSchOld:maxRayIndSchOld-1+N+N+fix(CP*N));%РЭЛЕЙ
% seq1Olddf(1,:) = seq1OldUpFcAwgnF0Dec500LpfFs1(1,401:401+N-1);%Теперь уже все по честному 
% seq12Olddf(1,:) = seq1OldUpFcAwgnF0Dec500LpfFs1(1,401:401-1+N+N+fix(CP*N));%АГБШ
% 
% seq1Olddf(2,:) = seq1OldUpFcRayF0Dec500LpfFs1(1,401:401-1+N);%Теперь уже все по честному
% seq12Olddf(2,:) = seq1OldUpFcRayF0Dec500LpfFs1(1,401:401-1+N+N+fix(CP*N));%РЭЛЕЙ
P1Old = sum(conj(seq1Olddf(1,1:N/2)).*seq1Olddf(1,1+N/2:N));
phiOld = angle(P1Old);
feAwgnOld = phiOld/pi;
CorrOld = exp(-1i*2*phiOld*t12/N);
seq12OlddfCorr = seq12Olddf(1,:).*CorrOld;
F1Old = fft(seq12OlddfCorr(1:N));F1Old = [F1Old(N/2+1:end),F1Old(1:N/2)];
F2Old = fft(seq12OlddfCorr(N+fix(CP*N)+1:end));F2Old = [F2Old(N/2+1:end),F2Old(1:N/2)];
F1F2Old = conj(F1Old).*F2Old;
% B(1) = power(abs(sum(F1F2.*conj(frs(1:2:end)))),2)/2/power(sum(power(F2,2)),2);
for i = 1:N/2
    F1F2shOld = circshift(F1F2Old,[0,2*(i-1)]);
    BgAwgnOld(i) = power(abs(sum(  F1F2shOld(1:2:end).*  conj(VkOld)      )),2)/2/power(sum(power(abs(F2Old),2)),2);
end
P1Old = sum(conj(seq1Olddf(2,1:N/2)).*seq1Olddf(2,1+N/2:N));
phiOld = angle(P1Old);
feRayOld = phiOld/pi;
CorrOld = exp(-1i*2*phiOld*t12/N);
seq12OlddfCorr = seq12Olddf(2,:).*CorrOld;
F1Old = fft(seq12OlddfCorr(1:N));F1Old = [F1Old(N/2+1:end),F1Old(1:N/2)];
F2Old = fft(seq12OlddfCorr(N+fix(CP*N)+1:end));F2Old = [F2Old(N/2+1:end),F2Old(1:N/2)];
F1F2Old = conj(F1Old).*F2Old;
% B(1) = power(abs(sum(F1F2.*conj(frs(1:2:end)))),2)/2/power(sum(power(F2,2)),2);
for i = 1:N/2
    F1F2shOld = circshift(F1F2Old,[0,2*(i-1)]);
    BgRayOld(i) = power(abs(sum(  F1F2shOld(1:2:end).*  conj(VkOld)      )),2)/2/power(sum(power(abs(F2Old),2)),2);
end
FdSchmidlAwgnOld = (1+256-find(BgAwgnOld==max(BgAwgnOld)))*2 +feAwgnOld;% 2*(513-(find(BgAwgn==max(BgAwgn)))) + feAwgn
FdSchmidlRayOld = (1+256-find(BgRayOld==max(BgRayOld)))*2 +feRayOld;%2*(513-(find(BgRay==max(BgRay)))) + feRay
if abs(FdSchmidlRayOld-FreqOffset )>MaxFreqError
    SchmidlOldError = 1;
    fprintf('ERROR, Schmidl Old!!!\n');
end

if EnableOutput
   fprintf('Error of proposed method is %.5f (awgn), %.5f (rayleygh)\n',FdProposedAwgn-FreqOffset,FdProposedRay-FreqOffset);
   fprintf('Error of modified Schmidl method is %.5f (awgn), %.5f (rayleygh)\n',FdSchmidlAwgn-FreqOffset,FdSchmidlRay-FreqOffset);      
   fprintf('Error of Schmidl method is %.5f (awgn), %.5f (rayleygh)\n',FdSchmidlAwgnOld-FreqOffset,FdSchmidlRayOld-FreqOffset);      
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
toc;
