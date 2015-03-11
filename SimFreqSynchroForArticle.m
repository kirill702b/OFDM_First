clear all; close all;


Ntests= 250;
SNR_n = 0;
FreqOffset_n = 0;
SNRs = 0:3:30;% 10.^((0:3:25)/10);
Freqs = 12.4;% 1.1:0.3:2.3;
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
            if ~ProposedError && ~ShmidlError
                FdSchmidlAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FdSchmidlAwgn;
                FdProposedAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FdProposedAwgn;
                FdSchmidlAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FdSchmidlRay;
                FdProposedAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FdProposedRay;            
            else
                FdSchmidlAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FreqOffset;
                FdProposedAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FreqOffset;
                FdSchmidlAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FreqOffset;
                FdProposedAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FreqOffset;
            end
            fprintf('SNR is %d, FreqOffset is %3.1f, test num is %d of %d; Proposed is better for %3.4f,%3.4f;\n',...
                SNR,FreqOffset,NumOfTest,Ntests,abs(FdSchmidlAwgn-FreqOffset)-abs(FdProposedAwgn-FreqOffset),abs(FdSchmidlRay-FreqOffset)-abs(FdProposedRay-FreqOffset));
                %             fprintf('SNR is %d, FreqOffset is %3.1f, test num is %d of %d, Errors:\n Proposed is %3.4f,%3.4f;\n  Schmidl is %3.4f,%3.4f.\n',...
%                 SNR,FreqOffset,NumOfTest,Ntests,FdProposedAwgn-FreqOffset,FdProposedRay-FreqOffset,FdSchmidlAwgn-FreqOffset,FdSchmidlRay-FreqOffset);
        end
        close all;
%         fprintf('Frequency offcet is %3.1f\n',FreqOffset);
    end  
    FreqOffset_n = 0;
    
%     fprintf('SNR is %d\n',SNR);
end
save('file_SNR_0_3_30_Ntests_500_n15_fd277div10.mat','FdSchmidlAll_awgn','FdProposedAll_awgn','FdSchmidlAll_ray','FdProposedAll_ray','Freqs','SNRs','Ntests');

