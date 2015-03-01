clear;
load('file_SNR_m30_30_30_Ntests_30.mat');

% FdSchmidlAll_awgn_mean =mean(FdSchmidlAll_awgn,3);
% FdProposedAll_awgn_mean = mean(FdProposedAll_awgn,3);
% 
% FdSchmidlAll_ray_mean =mean(FdSchmidlAll_ray,3);
% FdProposedAll_ray_mean = mean(FdProposedAll_ray,3);

% 
% figure;hold on;grid on; 
% plot(Freqs,FdSchmidlAll_awgn_mean(1,:));
% plot(Freqs,FdProposedAll_awgn_mean(1,:),'r');
% plot(Freqs,FdSchmidlAll_ray_mean(1,:),'g'); 
% plot(Freqs,FdProposedAll_ray_mean(1,:),'k');
% plot(Freqs,Freqs,'m');
% legend('FdSchmidlAll awgn','FdProposedAll awgn','FdSchmidlAll ray','FdProposedAll ray');

FdSchmidlAll_awgnSnr(:,:) = FdSchmidlAll_awgn(:,1,:);
FdProposedAll_awgnSnr(:,:) = FdProposedAll_awgn(:,1,:);

FdSchmidlAll_raySnr(:,:)=FdSchmidlAll_ray(:,1,:);
FdProposedAll_raySnr(:,:) = FdProposedAll_ray(:,1,:);

FdSchmidlAll_awgnSnr1 = sqrt(mean(power((FdSchmidlAll_awgnSnr-Freqs),2),2));
FdProposedAll_awgnSnr1 = sqrt(mean(power((FdProposedAll_awgnSnr-Freqs),2),2));
FdSchmidlAll_raySnr1 = sqrt(mean(power((FdSchmidlAll_raySnr-Freqs),2),2));
FdProposedAll_raySnr1 = sqrt(mean(power((FdProposedAll_raySnr-Freqs),2),2));

figure(1);hold on;grid on;
plot(SNRs,FdSchmidlAll_awgnSnr1);
plot(SNRs,FdProposedAll_awgnSnr1,'r');

figure(2);hold on;grid on;
plot(SNRs,FdSchmidlAll_raySnr1);
plot(SNRs,FdProposedAll_raySnr1,'r');