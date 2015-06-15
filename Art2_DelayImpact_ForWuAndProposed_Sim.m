clear all; close all;
FreqDop = 0.000;
Ntests = 1000;
SNR_n = 0;
FreqOffset_n = 0;
SNRs =10;%50;% 10.^((0:3:25)/10);
Freqs = 12.4;% 1.1:0.3:2.3;
delNs = 0:10:100;%150;
delN_n = 0;
dev = 50;
Nfft=2^17;
N=1024;%number of bits in preamble
Nnoise = 3;
frs3 = (randi(2,[1 N])-1.5)*2;% вторая последовательность
BP = ifft(frs3,N);
maxAwgnIndPro = Nnoise+1;
FdWuAwgn = zeros(length(SNRs),length(Freqs),length(delNs),Ntests);
FdMyAwgn = zeros(length(SNRs),length(Freqs),length(delNs),Ntests);
FdWuRay = zeros(length(SNRs),length(Freqs),length(delNs),Ntests);
FdMyRay = zeros(length(SNRs),length(Freqs),length(delNs),Ntests);
tic;
for SNR_n=1:length(SNRs)%SNR = SNRs    
    for FreqOffset_n = 1:length(Freqs)%FreqOffset = Freqs
        for delN_n = 1:length(delNs)
            fprintf('delN is %d\n',delNs(delN_n));
            parfor NumOfTest = 1:Ntests
                [seqAwgn,seqRay]=SequenceThroughChannelFunction(BP,N,Freqs(FreqOffset_n),SNRs(SNR_n),Nnoise,FreqDop,delNs(delN_n));
                FdWuAwgn(SNR_n,FreqOffset_n,delN_n,NumOfTest) = WuMethod(seqAwgn(1,maxAwgnIndPro:maxAwgnIndPro-1+N),BP);
                FdMyAwgn(SNR_n,FreqOffset_n,delN_n,NumOfTest) = MyMethod(seqAwgn(1,maxAwgnIndPro:maxAwgnIndPro-1+N),BP,Nfft,dev);
                FdWuRay(SNR_n,FreqOffset_n,delN_n,NumOfTest) = WuMethod(seqRay(1,maxAwgnIndPro:maxAwgnIndPro-1+N),BP);
                FdMyRay(SNR_n,FreqOffset_n,delN_n,NumOfTest) = MyMethod(seqRay(1,maxAwgnIndPro:maxAwgnIndPro-1+N),BP,Nfft,dev);
            end            
        end
    end    
    fprintf('SNR is %d\n',SNRs(SNR_n));
end
toc;
%%
save(['file_SNR_',num2str(SNRs(1)),'_',num2str(SNRs(end)),'_delNs_',num2str(delNs(1)),'_',num2str(delNs(end)),'_Ntests_',num2str(Ntests),'_Fdop_',num2str(FreqDop),'Hz','.mat'],'FdWuAwgn','FdMyAwgn','FdWuRay','FdMyRay','Freqs','SNRs','Ntests','delNs');

