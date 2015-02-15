clear all; close all;


Ntests= 10;
SNR_n = 0;
FreqOffset_n = 0;
SNRs = 10:10:40;% 10.^((0:3:25)/10);
Freqs = 12.4;% 1.1:0.3:2.3;
FdSchmidlAll = zeros(length(SNRs),length(Freqs),Ntests);
FdProposedAll = zeros(length(SNRs),length(Freqs),Ntests);
for SNR = SNRs
    SNR_n = SNR_n + 1;
    for FreqOffset = Freqs
        FreqOffset_n = FreqOffset_n + 1;
        for NumOfTest = 1:Ntests
            FreqSynchroForArticle
            FdSchmidlAll(SNR_n,FreqOffset_n,NumOfTest)=FdSchmidl;
            FdProposedAll(SNR_n,FreqOffset_n,NumOfTest)=FdProposed;
        end
%         FreqOffset
    end  
    FreqOffset_n = 0;
    SNR
end
%%

FdSchmidlAll1 =mean(FdSchmidlAll,3);
FdProposedAll1 = mean(FdProposedAll,3);

FdSchmidlErr = FdSchmidlAll-12.4;
FdSchmidlMSE = sqrt(mean(FdSchmidlErr.*FdSchmidlErr,3));

FdProposedErr = FdProposedAll-12.4;
FdProposedMSE = sqrt(mean(FdProposedErr.*FdProposedErr,3));


% figure;stem(Freqs,FdSchmidlAll1(1,:));grid on; 
% figure;stem(Freqs,FdProposedAll1(1,:));grid on;
