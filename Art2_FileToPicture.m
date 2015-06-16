clear; close all;
load('file_SNR_-10_-10_delNs_0_100_Ntests_1000_Fdop_0Hz.mat');

ErrWuAwgn = zeros(length(SNRs),length(delNs));
ErrMyAwgn = zeros(length(SNRs),length(delNs));
ErrWuRay = zeros(length(SNRs),length(delNs));
ErrMyRay = zeros(length(SNRs),length(delNs));


for i= 1:length(SNRs)
    for j = 1:length(delNs)
        ind = find(abs(FdWuAwgn((i),1,(j),:)-Freqs) >=0.5 );
        ErrWuAwgn(i,j)  = length(ind)/length(FdWuAwgn((i),1,(j),:));
        ind = find(abs(FdMyAwgn((i),1,(j),:)-Freqs) >=0.5 );
        ErrMyAwgn(i,j)  =length(ind)/length(FdWuAwgn((i),1,(j),:));
        ind = find(abs(FdWuRay((i),1,(j),:)-Freqs) >=0.5 );
        ErrWuRay(i,j)  = length(ind)/length(FdWuAwgn((i),1,(j),:));
        ind = find(abs(FdMyRay((i),1,(j),:)-Freqs) >=0.5 );
        ErrMyRay(i,j)  = length(ind)/length(FdWuAwgn((i),1,(j),:));
%         ErrWuAwgn(i,j)  = sqrt(mean(power(FdWuAwgn((i),1,(j),:)-Freqs,2)));
%         ErrMyAwgn(i,j)  = sqrt(mean(power(FdMyAwgn((i),1,(j),:)-Freqs,2)));
%         ErrWuRay(i,j)  = sqrt(mean(power(FdWuRay((i),1,(j),:)-Freqs,2)));
%         ErrMyRay(i,j)  = sqrt(mean(power(FdMyRay((i),1,(j),:)-Freqs,2)));
    end
    
    figure(i);
    plot(delNs/100,20*log10(ErrWuAwgn(i,:)));hold on;
    plot(delNs/100,20*log10(ErrMyAwgn(i,:)),'r');
    plot(delNs/100,20*log10(ErrWuRay(i,:)),'g');
    plot(delNs/100,20*log10(ErrMyRay(i,:)),'k');
%     plot(delNs/100,ErrWuAwgn(i,:));hold on;
%     plot(delNs/100,ErrMyAwgn(i,:),'r');
%     plot(delNs/100,ErrWuRay(i,:),'g');
%     plot(delNs/100,ErrMyRay(i,:),'k');
%     semilogy(delNs/100,ErrWuAwgn(i,:));hold on;
%     semilogy(delNs/100,ErrMyAwgn(i,:),'r');
%     semilogy(delNs/100,ErrWuRay(i,:),'g');
%     semilogy(delNs/100,ErrMyRay(i,:),'k');
    grid on;xlabel('Delay, ind');ylabel('MSR error, ici');legend('WuAwgn','MyAwgn','WuRay','MyRay');
end


