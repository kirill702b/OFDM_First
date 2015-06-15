clear; close all;
SNR = -10
dev = 50;
Nfft=2^17;
FreqOffset = 21.4;
N=1024;%number of bits in preamble
EnableGraphs = 0;
CP = 0.1; % Cyclic prefix length as part of symbol period
Nnoise = 3;
FreqDop = 0;
delN = 0;
frs3 = (randi(2,[1 N])-1.5)*2;% געמנאÿ ןמסכוהמגאעוכüםמסעü
BP = ifft(frs3,N);
tic;
[seqAwgn,seqRay]=SequenceThroughChannelFunction(BP,N,FreqOffset,SNR,Nnoise,FreqDop,delN);
toc;
if EnableGraphs
    figure;plot(abs(seqAwgn));hold on;plot(abs(seqRay),'m');plot(Nnoise+(1:length(BP)),abs(BP),'r');
end
%% Proposed
maxAwgnIndPro = Nnoise+1;
seq12df(1,:) = seqAwgn(1,maxAwgnIndPro:maxAwgnIndPro-1+N);%ÀÃÁØ
FdWuAwgn = WuMethod(seq12df(1,1:N),BP)
FdMyAwgn = MyMethod(seq12df(1,1:N),BP,Nfft,dev)

maxAwgnIndPro = Nnoise+1;
seq12df(2,:) = seqRay(1,maxAwgnIndPro:maxAwgnIndPro-1+N);%ÀÃÁØ
FdWuRay = WuMethod(seq12df(2,1:N),BP)
FdMyRay = MyMethod(seq12df(2,1:N),BP,Nfft,dev)
