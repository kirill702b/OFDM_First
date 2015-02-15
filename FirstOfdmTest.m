clear all; close all;
Nc = 1000; % Number of subcarriers
N = 1024; % Number of IFFT points
Rd = 18e6; % Data rate bit per second
CP = 0.1; % Cyclic prefix length as part of symbol period
FreqOffset = 12.4; % Frequency offcet in subcarrier spacing
NumOfOsc = 16; % Number of Oscillator for Jakes' model
Fc = 2e9; % Hz, Carrier frequency
OsFac = 10; %oversampling factor
Fs = Fc*OsFac; %sampling frequency
VehSpeed = 150; % kmph, Speed of vehicle
NumOfPath = 4; % Number of paths
DelOfPath = [0, 1, 2, 3]; % Delays of the paths in samples
AvRecPwr = [0, -3, -6, -9]; % Average received power, dB
Nnoise = 3000;%number of noise samples
NumOfMeth = 4;
phiRand = 2*pi/3;
%BPSK
%Training sequence


TrainSeqF1 = mseq(2,log2(N/2));
TrainSeqF1= [TrainSeqF1(1:Nc/2/2);zeros(N/2-Nc/2,1);TrainSeqF1(Nc/2/2+1:Nc/2) ];
TrainSeqT1 = ifft(TrainSeqF1,N/2).';

TrainSeqF2 = mseq(2,log2(N/2/2));
TrainSeqF2= [TrainSeqF2(1:Nc/2/2/2);zeros(N/2/2-Nc/2/2,1);TrainSeqF2(Nc/2/2/2+1:Nc/2/2) ];
TrainSeqT2 = ifft(TrainSeqF2,N/2/2).';



% TrainSeqF = mseq(2,log2(N/2));
% TrainSeqF = [TrainSeqF(1:Nc/2/2);zeros(N/2-Nc/2,1);TrainSeqF(Nc/2/2+1:Nc/2) ];
% TrainSeqT = ifft(TrainSeqF,N/2).';
% 
% TrainSeqF2 = mseq(2,log2(N/4));
% TrainSeqF2 = [TrainSeqF(1:Nc/2/2/2);zeros(N/2/2-Nc/2/2,1);TrainSeqF(Nc/2/2/2+1:Nc/2/2) ];
% TrainSeqT2 = ifft(TrainSeqF,N/2/2).';

TrainSeqT(1,1:N) = [TrainSeqT1(1:N/2),TrainSeqT1(1:N/2)];%Schmidl
TrainSeqT(2,1:N) = [TrainSeqT2(1:N/4),TrainSeqT2(1:N/4),-TrainSeqT2(1:N/4),-TrainSeqT2(1:N/4)];%Minn
TrainSeqT(3,1:N) = [TrainSeqT2(1:N/4),TrainSeqT2(N/4:-1:1),conj(TrainSeqT2(1:N/4)),conj(TrainSeqT2(N/4:-1:1))];%Park

% TrainSeqT(1,1:N) = [TrainSeqF(1:N/2)',TrainSeqF(1:N/2)'];%Schmidl
% TrainSeqT(2,1:N) = [TrainSeqF(1:N/4)',TrainSeqF(1:N/4)',-TrainSeqF(1:N/4)',-TrainSeqF(1:N/4)'];%Minn
% TrainSeqT(3,1:N) = [TrainSeqF(1:N/4)',TrainSeqF(N/4:-1:1)',conj(TrainSeqF(1:N/4))',conj(TrainSeqF(N/4:-1:1))'];%Park


% TrainSeqT(2,1:N) = [TrainSeqT2(1:N/4),TrainSeqT2(1:N/4),-TrainSeqT2(1:N/4),-TrainSeqT2(1:N/4)];%Minn
% TrainSeqT(3,1:N) = [TrainSeqT2(1:N/4),TrainSeqT2(N/4:-1:1),conj(TrainSeqT2(1:N/4)),conj(TrainSeqT2(N/4:-1:1))];%Park
TrainSeqT(4,1:N) = [zeros(1,N)];%Seung

RealSequence(1:NumOfMeth,1:2*Nnoise+N) = [zeros(NumOfMeth,Nnoise),TrainSeqT(1:NumOfMeth,1:N),zeros(NumOfMeth,Nnoise)];
% RealSequence(1:NumOfMeth,1:2*Nnoise+N) = [randn(NumOfMeth,Nnoise)+1i*randn(NumOfMeth,Nnoise),TrainSeqT(1:NumOfMeth,1:N),randn(NumOfMeth,Nnoise)+1i*randn(NumOfMeth,Nnoise)];
% RealSequence(1:NumOfMeth,1:2*Nnoise+N) = [0.001*ones(NumOfMeth,Nnoise),TrainSeqT(1:NumOfMeth,1:N),0.001*ones(NumOfMeth,Nnoise)];
RealSequenceNoised = awgn(RealSequence,50);


t = (1:(2*Nnoise+N))/Fs;
for i=1:NumOfMeth
    RealSignal(i,1:2*Nnoise+N) = RealSequence(i,1:2*Nnoise+N).*exp(1i*2*pi*Fc*t);
end




RayleighCh = comm.RayleighChannel('PathDelays',DelOfPath,'AveragePathGains',AvRecPwr);
for i=1:NumOfMeth
    RealSignalCh(i,:) = step(RayleighCh,awgn(RealSignal(i,:),50).');
    reset(RayleighCh);
end

t1 = (1:length(RealSignalCh(1,:)))/Fs;
for i=1:NumOfMeth
    SeqRec(i,:) = RealSignalCh(i,:).*exp(-1i*(2*pi*Fc*t1+phiRand));
end

% RealSequenceAwgn(:,:) = awgn(RealSequence(:,:),1000);
% for j=1:length(RealSequence(1,:))-N
%     P1 = sum(conj(RealSequenceAwgn(1,j:j+N/2-1)).*RealSequenceAwgn(1,j+N/2:j+N/2-1+N/2));
%     R1 = sum(power(abs(RealSequenceAwgn(1,j+N/2:j+N/2-1+N/2)),2));
%     RespOfFindIn(1,j) = power(abs(P1),2)/power(R1,2);
% end
% 
% L2 = N/4;
% for j=1:length(RealSequenceAwgn(2,:))-N
%     P1 = sum(conj(RealSequenceAwgn(2,j:j+L2-1)).*RealSequenceAwgn(2,j+L2:j+L2-1+L2)) + sum(conj(RealSequenceAwgn(2,j+2*L2:j+L2-1+2*L2)).*RealSequenceAwgn(2,j+L2+2*L2:j+L2-1+L2+2*L2));
%     R1 = sum(power(abs(RealSequenceAwgn(2,j+L2:j+L2-1+L2)),2)) + sum(power(abs(RealSequenceAwgn(2,j+L2+2*L2:j+L2-1+L2+2*L2)),2));
%     RespOfFindIn(2,j) = power(abs(P1),2)/power(R1,2);
% end



figure;hold on;plot(real(RealSequenceNoised(1,:)));
% plot(abs(RealSequenceAwgn(1,:)),'r');

L1= N/2;
for j=1:length(SeqRec(1,:))-N
    P1 = sum(conj(SeqRec(1,j:j+L1-1)).*SeqRec(1,j+L1:j+L1-1+L1));
    R1 = sum(power(abs(SeqRec(1,j+L1:j+L1-1+L1)),2));
    RespOfFind(1,j) = power(abs(P1),2)/power(R1,2);
end
L2 = N/4;
for j=1:length(SeqRec(2,:))-N
    P1 = sum(conj(SeqRec(2,j:j+L2-1)).*SeqRec(2,j+L2:j+L2-1+L2)) + sum(conj(SeqRec(2,j+2*L2:j+L2-1+2*L2)).*SeqRec(2,j+L2+2*L2:j+L2-1+L2+2*L2));
    R1 = sum(power(abs(SeqRec(2,j+L2:j+L2-1+L2)),2)) + sum(power(abs(SeqRec(2,j+L2+2*L2:j+L2-1+L2+2*L2)),2));
    RespOfFind(2,j) = power(abs(P1),2)/power(R1,2);
end

%% Tests
for j=1:length(RealSequenceNoised(1,:))-N
    P1 = sum(conj(RealSequenceNoised(1,j:j+N/2-1)).*RealSequenceNoised(1,j+N/2:j+N/2-1+N/2));
    R1 = sum(power(abs(RealSequenceNoised(1,j+N/2:j+N/2-1+N/2)),2));
    RespOfFindIn(1,j) = power(abs(P1),2)/power(R1,2);
end

L2 = N/4;
for j=1:length(RealSequenceNoised(2,:))-N
    P1 = sum(conj(RealSequenceNoised(2,j:j+L2-1)).*RealSequenceNoised(2,j+L2:j+L2-1+L2)) ...
        + sum(conj(RealSequenceNoised(2,j+2*L2:j+L2-1+2*L2)).*RealSequenceNoised(2,j+L2+2*L2:j+L2-1+L2+2*L2));
    R1 = sum(power(abs(RealSequenceNoised(2,j+L2:j+L2-1+L2)),2)) + sum(power(abs(RealSequenceNoised(2,j+L2+2*L2:j+L2-1+L2+2*L2)),2));
    RespOfFindIn(2,j) = power(abs(P1),2)/power(R1,2);
end
L3 = N/2;
for j=N/2+1:length(RealSequenceNoised(2,:))-N/2
    P1 = 0;R1 = 0;
    for k = 0:L3
        P1 = P1 + RealSequenceNoised(2,j+k)*RealSequenceNoised(2,j-k);
        R1 = R1 + power(abs(RealSequenceNoised(2,j+k)),2);
    end  
    RespOfFindIn(3,j-N/2) = power(abs(P1),2)/power(R1,2);
end
%%

figure;hold on;
plot(RespOfFindIn(1,:));
plot(RespOfFindIn(2,:),'r');
plot(RespOfFindIn(3,:),'g');
% 
% figure;hold on;
% plot(RespOfFind(1,:));
% plot(RespOfFind(2,:),'r');







