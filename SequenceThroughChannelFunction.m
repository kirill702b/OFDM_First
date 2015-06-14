function [ seqInUpFcChF0PreDecLpfFs1 ] = SequenceThroughChannelFunction( Preamble,N,FreqOffset,SNR,isRayOrAwgnOnly,Nnoise,FreqDop,delN)
    Rd = 18e6; % Data rate bit per second
    CP = 0.1; % Cyclic prefix length as part of symbol period
    Fs1 =(Rd/(1-CP));%20mHz, 50ns 
    Fc = 2e8; % Hz, Carrier frequency 2GHz
    OsFac = 10; %oversampling factor
    Fs = Fc*OsFac; %sampling frequency Fs = 20GHz
    NumOfPath = 4; % Number of paths
    DelOfPath = [0, 1, 2, 3]/Fs; % Delays of the paths in samples
    AvRecPwr = [0, -3, -6, -9]; % Average received power, dB
    seqIn = [zeros(1,Nnoise),real(Preamble),zeros(1,Nnoise)]+1i* [zeros(1,Nnoise),imag(Preamble),zeros(1,Nnoise)];
    seqInUp = resample(seqIn,Fs/Fs1,1); 
    tUp3 = (1:length(seqInUp))/Fs;
    seqInUpFc = seqInUp.*exp(1i*2*pi*(Fc+FreqOffset*Fs1/N)*tUp3);
    OrderOfFil1 = 64;
    OrderOfFil2 = 128;
    if ~exist('RayCh1')
        RayCh1 = rayleighchan(1/Fs,FreqDop,DelOfPath,AvRecPwr);
        RayCh1.ResetBeforeFiltering = 1;
    end
    if isRayOrAwgnOnly
        seqInUpFcCh = awgn(filter(RayCh1,seqInUpFc),SNR);
    else
        seqInUpFcCh = awgn(seqInUpFc,SNR);%сигнал после АГБШ канала
    end
    seqInUpFcCh = [seqInUpFcCh(delN+1:end),zeros(1,delN)];
    seqInUpFcChF0 = seqInUpFcCh.*exp(-1i*2*pi*Fc*tUp3);
    DecFactor = Fs/Fs1;
    if DecFactor-fix(DecFactor)~=0
        display('Decimation factor is not integer!!! Exit.');
        return;
    end
    Decs = sort(factor(DecFactor),'descend');
    seqInUpFcChF0PreDec = decimate(seqInUpFcChF0,Decs(1),OrderOfFil1,'fir');
    for j=2:length(Decs)-1
        seqInUpFcChF0PreDec = decimate(seqInUpFcChF0PreDec,Decs(j),OrderOfFil1,'fir');
    end
    seqInUpFcChF0PreDecLpf=seqInUpFcChF0PreDec;
    seqInUpFcChF0PreDecLpfFs1 = decimate(seqInUpFcChF0PreDecLpf,2,OrderOfFil2,'fir');
end

