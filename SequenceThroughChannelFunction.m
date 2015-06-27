function [ seqOutAwgn, seqOutRay ] = SequenceThroughChannelFunction( Preamble,N,FreqOffset,SNR,Nnoise,FreqDop,delN)
    HaveCarrier = 1;
    Rd = 18e6; % Data rate bit per second
    CP = 0.1; % Cyclic prefix length as part of symbol period
    Fs1 =(Rd/(1-CP));%20mHz, 50ns 
    Fc = 2e8; % Hz, Carrier frequency 2GHz
    OsFac = 10; %oversampling factor
    Fs = Fc*OsFac; %sampling frequency Fs = 20GHz
    NumOfPath = 4; % Number of paths    
    DelOfPath = [0, 1, 2, 3]/Fs;    
    AvRecPwr = [0, -3, -6, -9]; % Average received power, dB
    seqIn = [zeros(1,Nnoise),real(Preamble),zeros(1,Nnoise)]+1i* [zeros(1,Nnoise),imag(Preamble),zeros(1,Nnoise)];
    if HaveCarrier
        seqInUp = resample(seqIn,Fs/Fs1,1); 
        tUp3 = (1:length(seqInUp))/Fs;
    else
        seqInUp = seqIn; 
        tUp3 = (1:length(seqInUp))/Fs1;
    end
    if HaveCarrier
        seqInUpFc = seqInUp.*exp(1i*2*pi*(Fc+FreqOffset*Fs1/N)*tUp3);
    else
        seqInUpFc = seqInUp.*exp(1i*2*pi*(FreqOffset*Fs1/N)*tUp3);
    end
    OrderOfFil1 = 64;
    OrderOfFil2 = 128;
    if HaveCarrier
        if ~exist('RayCh1')
            RayCh1 = rayleighchan(1/Fs,FreqDop,DelOfPath,AvRecPwr);
            RayCh1.ResetBeforeFiltering = 1;
        end
    else
        if ~exist('RayCh1')
            RayCh1 = rayleighchan(1/Fs1,FreqDop,DelOfPath,AvRecPwr);
            RayCh1.ResetBeforeFiltering = 1;
        end
    end
    if SNR~=100
        seqInUpFcAwgn = awgn(seqInUpFc,SNR,'measured');
    else
        seqInUpFcAwgn =seqInUpFc;
    end
    seqInUpFcAwgn = [seqInUpFcAwgn(delN+1:end),zeros(1,delN)];
    if HaveCarrier
        seqInUpFcAwgnF0 = seqInUpFcAwgn.*exp(-1i*2*pi*Fc*tUp3);
    else
        seqInUpFcAwgnF0 = seqInUpFcAwgn;
    end
    DecFactor = Fs/Fs1;
    if DecFactor-fix(DecFactor)~=0
        display('Decimation factor is not integer!!! Exit.');
        return;
    end
    Decs = sort(factor(DecFactor),'descend');
    if HaveCarrier
        seqInUpFcAwgnF0PreDec = decimate(seqInUpFcAwgnF0,Decs(1),OrderOfFil1,'fir');
        for j=2:length(Decs)-1
            seqInUpFcAwgnF0PreDec = decimate(seqInUpFcAwgnF0PreDec,Decs(j),OrderOfFil1,'fir');
        end
        seqInUpFcAwgnF0PreDecLpf=seqInUpFcAwgnF0PreDec;   
        seqOutAwgn = decimate(seqInUpFcAwgnF0PreDecLpf,2,OrderOfFil2,'fir');
    else
        seqInUpFcAwgnF0PreDec = seqInUpFcAwgnF0;
        for j=2:length(Decs)-1
            seqInUpFcAwgnF0PreDec = seqInUpFcAwgnF0PreDec;
        end
        seqInUpFcAwgnF0PreDecLpf=seqInUpFcAwgnF0PreDec;   
        seqOutAwgn = seqInUpFcAwgnF0PreDecLpf;
    end
    
    
    if SNR~=100
        seqInUpFcRay =awgn(filter(RayCh1,seqInUpFc),SNR,'measured') ;%сигнал после АГБШ канала
    else
        seqInUpFcRay=filter(RayCh1,seqInUpFc);
    end
    seqInUpFcRay = [seqInUpFcRay(delN+1:end),zeros(1,delN)];
    if HaveCarrier
        seqInUpFcRayF0 = seqInUpFcRay.*exp(-1i*2*pi*Fc*tUp3);
    else
        seqInUpFcRayF0 = seqInUpFcRay;
    end
    DecFactor = Fs/Fs1;
    if DecFactor-fix(DecFactor)~=0
        display('Decimation factor is not integer!!! Exit.');
        return;
    end
    Decs = sort(factor(DecFactor),'descend');
    if HaveCarrier
        seqInUpFcRayF0PreDec = decimate(seqInUpFcRayF0,Decs(1),OrderOfFil1,'fir');
        for j=2:length(Decs)-1
            seqInUpFcRayF0PreDec = decimate(seqInUpFcRayF0PreDec,Decs(j),OrderOfFil1,'fir');
        end
        seqInUpFcRayF0PreDecLpf=seqInUpFcRayF0PreDec;   
        seqOutRay = decimate(seqInUpFcRayF0PreDecLpf,2,OrderOfFil2,'fir');
    else
        seqInUpFcRayF0PreDec = seqInUpFcRayF0;
        for j=2:length(Decs)-1
            seqInUpFcRayF0PreDec = seqInUpFcRayF0PreDec;
        end
        seqInUpFcRayF0PreDecLpf=seqInUpFcRayF0PreDec;   
        seqOutRay = seqInUpFcRayF0PreDecLpf;
    end
end

