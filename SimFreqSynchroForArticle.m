clear all; close all;


Ntests= 20;
SNR_n = 0;
FreqOffset_n = 0;
SNRs = 15:5:30;% 10.^((0:3:25)/10);
Freqs = 8.3;% 1.1:0.3:2.3;
FdSchmidlAll_awgn = zeros(length(SNRs),length(Freqs),Ntests);
FdProposedAll_awgn = zeros(length(SNRs),length(Freqs),Ntests);
FdSchmidlAll_ray = zeros(length(SNRs),length(Freqs),Ntests);
FdProposedAll_ray = zeros(length(SNRs),length(Freqs),Ntests);
for SNR = SNRs
    SNR_n = SNR_n + 1;
    for FreqOffset = Freqs
        FreqOffset_n = FreqOffset_n + 1;
        for NumOfTest = 1:Ntests
            FreqSynchroForArticle_withRF
            FdSchmidlAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FdSchmidlAwgn;
            FdProposedAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FdProposedAwgn;
            FdSchmidlAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FdSchmidlRay;
            FdProposedAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FdProposedRay;
            fprintf('Number of test is %d of %d\n',NumOfTest,Ntests);
        end
        close all;
        fprintf('Frequency offcet is %3.1f\n',FreqOffset);
    end  
    FreqOffset_n = 0;
    
    fprintf('SNR is %d\n',SNR);
end
save('file2.mat','FdSchmidlAll_awgn','FdProposedAll_awgn','FdSchmidlAll_ray','FdProposedAll_ray','Freqs','SNRs','Ntests');

