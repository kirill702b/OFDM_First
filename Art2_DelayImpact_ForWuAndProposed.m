SNR = 50;
dev = 50;
Nfft=2^17;
FreqOffset = 20.4;
N=1024;%number of bits in preamble
EnableGraphs = 0;
CP = 0.1; % Cyclic prefix length as part of symbol period
Nnoise = 400;
FreqDop = 0;
delN = 0;
frs3 = (randi(2,[1 N])-1.5)*2;% вторая последовательность
frs4 = (randi(2,[1 N/4])-1.5)*2;
C = ifft(frs4,N/4)/2;
BP = ifft(frs3,N);
CiCcCciC=[C,C(end:-1:1),conj(C),conj(C(end:-1:1))];
CiCcCciC_B = [CiCcCciC,zeros(1,fix(CP*N)),BP];%A,A];
GciG=CiCcCciC_B;%[CiCcCciC,zeros(1,fix(CP*N)),CiCcCciC];
G3=GciG(end-N+1:end);
tic;
seq3UpFcAwgnF0Dec500LpfFs1=SequenceThroughChannelFunction(CiCcCciC_B,N,FreqOffset,SNR,0,Nnoise,FreqDop,delN);
toc;
if EnableGraphs
    figure;plot(abs(seq3UpFcAwgnF0Dec500LpfFs1));hold on;plot(Nnoise+(1:length(CiCcCciC_B)),abs(CiCcCciC_B),'r');
end
%% Proposed
maxAwgnIndPro = Nnoise+1;%find(RespOfFind(3,:)==max(RespOfFind(3,:)));

seq12df(1,:) = seq3UpFcAwgnF0Dec500LpfFs1(1,maxAwgnIndPro:maxAwgnIndPro-1+N+N+fix(CP*N));%АГБШ
FdWu = WuMethod(seq12df(1,1+N+fix(CP*N):N+N+fix(CP*N)),G3)
FdMy = MyMethod(seq12df(1,1+N+fix(CP*N):N+N+fix(CP*N)),G3,Nfft,dev)
