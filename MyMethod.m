function [ FdOut ] = MyMethod( inSeq,RealSeq,Nfft,dev )
    N = length(inSeq);
    Fk = fft(conj(RealSeq).*inSeq,N);
    Fk = [Fk(N/2+1:end),Fk(1:N/2)];
    kmax =  find(abs(Fk)==max(abs(Fk)));
    Fcoarse = kmax-N/2;
    if kmax>=2 && kmax <=1023
        if abs(Fk(kmax-1))<= abs(Fk(kmax+1))
            alp = 1;
        else
            alp = -1;
        end    
    end
    fmax1 = ((Fcoarse-1-alp)/N*Nfft);
    fmax2 =  ((Fcoarse-1+alp)/N*Nfft);
    kk = 1;
    ts = conj(RealSeq).*inSeq.*window(@hamming,N).';%hamming
    for jj = min(fmax1,fmax2)-dev:max(fmax1,fmax2)+dev
        FHRAll(kk) = sum(ts.*exp(-1i*2*pi/Nfft*jj*(0:N-1)));
        kk=kk+1;
    end
    FHR = FHRAll(dev+1:end-dev);
    indm = find(abs(FHR(1:end)) == max(abs(FHR)))+dev;
    indmm = sum(sqrt(abs(FHRAll(indm-dev:indm+dev))).*(indm-dev:indm+dev)) / sum(sqrt(abs(FHRAll(indm-dev:indm+dev))));
    FdOut = Fcoarse-2 + (indmm-1 -dev)*N/Nfft;
end

