clear all; close all;
FreqDop = 10;
Ntests=1000;
SNR_n = 0;
FreqOffset_n = 0;
SNRs = 0:3:36;% 10.^((0:3:25)/10);
Freqs = 12.4;% 1.1:0.3:2.3;
FdSchmidlAll_awgn = zeros(length(SNRs),length(Freqs),Ntests);
FdProposedAll_awgn = zeros(length(SNRs),length(Freqs),Ntests);
FdSchmidlAll_ray = zeros(length(SNRs),length(Freqs),Ntests);
FdProposedAll_ray = zeros(length(SNRs),length(Freqs),Ntests);
FdSchmidlOldAll_awgn = zeros(length(SNRs),length(Freqs),Ntests);
FdSchmidlOldAll_ray = zeros(length(SNRs),length(Freqs),Ntests);
for SNR = SNRs
    SNR_n = SNR_n + 1;
    for FreqOffset = Freqs
        FreqOffset_n = FreqOffset_n + 1;
        for NumOfTest = 1:Ntests
            FreqSynchroForArticle_withRF
%             if ~ProposedError && ~SchmidlError && ~SchmidlOldError
                FdSchmidlAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FdSchmidlAwgn;
                FdProposedAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FdProposedAwgn;
                FdSchmidlAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FdSchmidlRay;
                FdProposedAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FdProposedRay;    
                FdSchmidlOldAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FdSchmidlAwgnOld;
                FdSchmidlOldAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FdSchmidlRayOld;
%             else
%                 FdSchmidlAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FreqOffset;
%                 FdProposedAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FreqOffset;
%                 FdSchmidlAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FreqOffset;
%                 FdProposedAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FreqOffset;
%                 FdSchmidlOldAll_awgn(SNR_n,FreqOffset_n,NumOfTest)=FreqOffset;
%                 FdSchmidlOldAll_ray(SNR_n,FreqOffset_n,NumOfTest)=FreqOffset;
%             end
            fprintf('SNR is %d, FreqOffset is %3.1f, test num is %d of %d;\n Proposed Err      %3.4f %3.4f\n Schmidl New Error %3.4f %3.4f\n Schmidl Old Error %3.4f %3.4f\n',...
                SNR,FreqOffset,NumOfTest,Ntests,(FdProposedAwgn-FreqOffset),(FdProposedRay-FreqOffset),(FdSchmidlAwgn-FreqOffset),(FdSchmidlRay-FreqOffset),(FdSchmidlAwgnOld-FreqOffset),(FdSchmidlRayOld-FreqOffset));
                %             fprintf('SNR is %d, FreqOffset is %3.1f, test num is %d of %d, Errors:\n Proposed is %3.4f,%3.4f;\n  Schmidl is %3.4f,%3.4f.\n',...
%                 SNR,FreqOffset,NumOfTest,Ntests,FdProposedAwgn-FreqOffset,FdProposedRay-FreqOffset,FdSchmidlAwgn-FreqOffset,FdSchmidlRay-FreqOffset);
        end
        close all;
%         fprintf('Frequency offcet is %3.1f\n',FreqOffset);
    end  
    FreqOffset_n = 0;
    
%     fprintf('SNR is %d\n',SNR);
end
%%
save(['file_SNR_',num2str(SNRs(1)),'_',num2str(SNRs(end)),'_Ntests_',num2str(Ntests),'_Fdop_',num2str(FreqDop),'Hz','.mat'],'FdSchmidlAll_awgn','FdProposedAll_awgn','FdSchmidlAll_ray','FdProposedAll_ray','Freqs','SNRs','Ntests','FdSchmidlOldAll_awgn','FdSchmidlOldAll_ray');

