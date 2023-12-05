function [ShiftedSBEs,shift]=allignSpikeShapes(convM,Method)
%[ShiftedSBEs,meanSBE,stdSBE,shift]=allignSpikeShapes(convM,Method)
%calculate the average SBE from a given set of bursts (already convoluted)
%convM - the selected bursts convoluted matrix.
%convM(bursts,neurons,data_bins);
%ShiftedSBEs - The shifted bursts for the case of one neurons.
%Method='Plot' for plotting

%21.2.10 - Mark (modified temporarily for spike sorting from averageSBE

if nargin<2
    Method='Non Plot';
end

[BurstNum, NeuronNum, BurstSize]=size(convM);
Bursts=zeros(BurstNum,BurstSize*(2*NeuronNum+1));

%build a matrix of 1-dimensional vectors of the bursts:
%put zeros between each neuron recording, and at bagining and ending of 1st and last neurons.
if NeuronNum==1
        Bursts=[zeros(BurstNum,BurstSize) squeeze(convM) zeros(BurstNum,BurstSize)];
else
       % [0 0 0 neuron1 0 0 0 neuron2 0 0 0...neuronN 0 0 0]
       Bursts=[zeros(BurstSize,BurstNum);reshape(permute(cat(3,convM,zeros(BurstNum,NeuronNum,BurstSize)),[3 2 1]),NeuronNum*(BurstSize*2),BurstNum)]';
end


%calculate the first approximation:
fprintf('Building initial approximated burst...');
Burst1=Bursts(1,:);
shift=zeros(1,BurstNum);    %vector of shifts of each burst from 1st burst.
Group=Bursts;
step=0;
while size(Group,1)>1,
    step=step+1;
    for i=1:2:size(Group,1),
        if i+1<=size(Group,1),
            Burst1=Group(i,:);
            Burst2=Group(i+1,:);
            nanVals=isnan(Burst1);Burst1(nanVals)=mean(Burst1(~nanVals));nanVals=isnan(Burst2);Burst2(nanVals)=mean(Burst2(~nanVals)); %Verify that there are no NaNs
            [XCF, Lags, Bounds] = crosscorr2(Burst1, Burst2, BurstSize-1);
            %[XCF,Lags]=xcorr(Burst1, Burst2,BurstSize,'coeff') ;
            [MaxCrossCorrValue,MaxValueLocation]=max(XCF);
            for j=1:2^(step-1),
                if (2^(step-1)*i+j)<=BurstNum,
                    shift(2^(step-1)*i+j)=shift(2^(step-1)*i+j)+Lags(MaxValueLocation);
                end
            end
            %Burst1=Burst1(BurstSize+1:end-BurstSize) + Burst2(BurstSize+Lags(MaxValueLocation)+1:end-BurstSize+Lags(MaxValueLocation));
            if Lags(MaxValueLocation)>0,
                Burst1(1:end-Lags(MaxValueLocation))=Burst1(1:end-Lags(MaxValueLocation))+Burst2(Lags(MaxValueLocation)+1:end);
            else
                Burst1(1-Lags(MaxValueLocation):end)=Burst1(1-Lags(MaxValueLocation):end)+Burst2(1:end+Lags(MaxValueLocation));
            end
            %replace burst in summed bursts:
            Group(i,:)=Burst1;
        end
    end
    Group(2:2:end,:,:)=[];
end
Burst1=squeeze(Group(1,:,:));



%iterate for higher approximations. this is because the first bursts were
%positioned in a more robust manner.
ChangedLocations=ones(1,BurstNum); %when all zeros - then bursts are not changed in relative postion.
iteration=1;
fprintf('\nIterating (percent of correct positioned bursts): ');
while (iteration<=40 & sum(ChangedLocations)>0)
    fprintf('%d , ' , 100-round(sum(ChangedLocations)/BurstNum*100));
    for i=1:BurstNum,
        Burst2=Bursts(i,:);
        %remove the burst that is to be rechecked:
        if shift(i)>=0,
            Burst1(1:end-shift(i))=Burst1(1:end-shift(i))-Burst2(shift(i)+1:end);
        else
            Burst1(1-shift(i):end)=Burst1(1-shift(i):end)-Burst2(1:end+shift(i));
        end
        nanVals=isnan(Burst1);Burst1(nanVals)=mean(Burst1(~nanVals));nanVals=isnan(Burst2);Burst2(nanVals)=mean(Burst2(~nanVals)); %Verify that there are no NaNs
        %cross-correlate the two vectors:
        [XCF, Lags, Bounds] = crosscorr2(Burst1, Burst2, BurstSize-1);
        %[XCF,Lags]=xcorr(Burst1, Burst2,BurstSize,'coeff') ;
        %find the location of the max value in burst:
        [MaxCrossCorrValue,MaxValueLocation]=max(XCF);
        if Lags(MaxValueLocation)==shift(i),
            ChangedLocations(i)=0;
        else,
            ChangedLocations(i)=1;
            shift(i)=Lags(MaxValueLocation);
        end
        %calculate the mean burst from tht two shifted vectors:
        if shift(i)>0,
            Burst1(1:end-shift(i))=Burst1(1:end-shift(i)) + Burst2(shift(i)+1:end);
        else
            Burst1(1-shift(i):end)=Burst1(1-shift(i):end) + Burst2(1:end+shift(i));
        end

        %build the new burst for next iteration:
        %        Burst1=[zeros(1,BurstSize),Burst1,zeros(1,BurstSize)];
    end
    iteration=iteration+1;
end
fprintf('\ndone!\n');

%plot the mean and std bursts:
[ShiftedSBEs]=StdSBE(convM,shift,Method);

function [ShiftedSBEs]=StdSBE(convM,shift,Method)
%function [meanSBE,stdSBE]=StdSBE(convM,shift);
%calculates the mean and std of burst according to shift as calculated in
%AverageSBE
%convM - the selected bursts convoluted matrix.
%convM(bursts,neurons,data_bins);
%shift - the shift in bins of each burst relative to the mean
%ShiftedSBEs - The shifted bursts for the case of one neurons.

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
        Bursts(i,:)=[zeros(1,BurstSize), reshape([squeeze(convM(i,:,:)) , zeros(NeuronNum,BurstSize)]',1,NeuronNum*(BurstSize*2))];
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

if NeuronNum==1
    ShiftedSBEs=ShiftedBursts(:,BurstSize+1:BurstSize+BurstSize);
else
    for j=1:NeuronNum
        ShiftedSBEs(:,j,:)=ShiftedBursts(:,BurstSize*(j*2-1)+1:BurstSize*(j*2-1)+BurstSize);
    end
end

%reshape the final burst back to matrix format:
for i=1:NeuronNum,
    meanSBE(i,:)=Burst1(BurstSize*(2*i-1)+1:BurstSize*(2*i-1)+BurstSize);
    stdSBE(i,:)=Std1(BurstSize*(2*i-1)+1:BurstSize*(2*i-1)+BurstSize);
end

if strcmp(Method,'Plot')
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
    
end

