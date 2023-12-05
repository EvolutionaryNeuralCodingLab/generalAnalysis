function [meanSBE,stdSBE,ShiftedBursts]=StdSBE(convM,shift);
%function [meanSBE,stdSBE]=StdSBE(convM,shift);
%calculates the mean and std of burst according to shift as calculated in
%AverageSBE
%convM - the selected bursts convoluted matrix.
%convM(bursts,neurons,data_bins);
%shift - the shift in bins of each burst relative to the mean
%The shifted bursts for the case of one neurons.

%04.5.09 - Mark: based on AverageSBE;

[BurstNum, NeuronNum, BurstSize]=size(convM);

%build a matrix of 1-dimensional vectors of the bursts:
%put zeros between each neuron recording, and at bagining and ending of 1st and last neurons.
if NeuronNum==1
    for i=1:BurstNum,
        Bursts(i,:)=[zeros(1,BurstSize) squeeze(convM(i,:,:))' zeros(1,BurstSize)];
    end  
else
    for i=1:BurstNum,
        Bursts(i,:)=[zeros(1,BurstSize), reshape([squeeze(convM(i,:,:),1) , zeros(NeuronNum,BurstSize)]',1,NeuronNum*(BurstSize*2))];
    end           % [0 0 0 neuron1 0 0 0 neuron2 0 0 0...neuronN 0 0 0]
end

%build a matrix of the shifted bursts:
ShiftedBursts=zeros(size(Bursts));
for i=1:BurstNum,
    ShiftedBursts(i,(BurstSize+2):(end-BurstSize-shift(i)+1))=Bursts(i,(BurstSize+shift(i)+1):(end-BurstSize));
end
Burst1=mean(ShiftedBursts);
Std1=std(ShiftedBursts);
%Burst1=Burst1/BurstNum;

%reshape the final burst back to matrix format:
for i=1:NeuronNum,
    meanSBE(i,:)=Burst1(BurstSize*(2*i-1)+1:BurstSize*(2*i-1)+BurstSize);
    stdSBE(i,:)=Std1(BurstSize*(2*i-1)+1:BurstSize*(2*i-1)+BurstSize);
end

%plot the mean burst:
figure;
for i=1:NeuronNum,
    hold on;
    plot(meanSBE(NeuronNum-i+1,:)+i-1);
end
hold off;

%plot mean with (+) and (-) std of burst:
figure;
ZeroVec=zeros(1,length(meanSBE));
P=2.0;
for i=1:NeuronNum,
    hold on;
    plot(1:length(meanSBE),ZeroVec+i-1,'k--');
    plot(meanSBE(NeuronNum-i+1,:)/P+stdSBE(NeuronNum-i+1,:)/P+i-1,'r');
    plot(meanSBE(NeuronNum-i+1,:)/P-stdSBE(NeuronNum-i+1,:)/P+i-1,'r');
    plot(meanSBE(NeuronNum-i+1,:)/P+i-1,'b');
end
hold off;

