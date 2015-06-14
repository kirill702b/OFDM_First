function [ FdOut ] = WuMethod(inSeq,RealSeq)
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
        Ffine = alp/(abs(Fk(kmax))/abs(Fk(kmax+alp))+1);
    else
        Ffine = 0;    
    end
    FdOut = Fcoarse + Ffine - 1;
end

