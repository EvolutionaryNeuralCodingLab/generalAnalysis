% [R,A]=orderClustering(M,varargin)
% Function purpose : Calculates a similarity matrix between different activity profiles
%
% Function recives :   M  - The activity matrix [Bind x Neurons x Activity trace length]
%                                               varagin (list of possible options)
%
% Function give back :  R 
% Recommended usage: [R]=orderClustering(M);
%
% Last updated : 21/09/18
function [R,TimeOrder,RS,A,AS]=orderClustering(M,varargin)

input='activity'; %options: 'activity'/'delays'
Threshold=0.1; %threshold for insignificant delays
weights='Binary';%options 'NonBinary'/'Binary'
method='FirstSpike';
averageOverNeurons=false; %if false, will compare sequences to the average propagation for each neuron
usePenalty=true;
nShuffles=10;
nearestNeigh=10;


%Collects all options
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

if strcmp(input,'activity')
    % intialize variables
    [NBs NN L]=size(M); %[NB number X Neuron number X Activity window
    
    if strcmp(method,'FirstSpike')
        M2=zeros(size(M));
        for i=1:size(M,1)
            tmp=findfirst(M(i,:,:), 3);
            p=find(tmp>0);
            for j=1:numel(p)
                M2(i,p(j),tmp(p(j)))=1;
            end
        end
        M=M2;
    end
    
    TimeOrder=zeros(NBs,NN);
    A=zeros(NBs,NN,NN);
    AS=zeros(NBs,NN,NN,nShuffles);
    
    %Build a binary delay matrix of all neuronal pears for every NB (network burst)
    Times=(1:L)';
    for i=1:NBs
        tmpM=squeeze(M(i,:,:));
        TimeOrder(i,:)=tmpM*Times./sum(tmpM,2); %calculates the average time delay for every neuron
        TimeOrderShuff=TimeOrder(i,:);
        pNonNaN=find(~isnan(TimeOrder(i,:)));
        %Build a delay matrix for the ith NB
        
        for l=1:nShuffles
            TimeOrderShuff(pNonNaN)=TimeOrder(i,pNonNaN(randperm(numel(pNonNaN))));
            for m=1:NN
                for n=1:NN
                    A(i,m,n)=TimeOrder(i,m)-TimeOrder(i,n);
                    AS(i,m,n,l)=TimeOrderShuff(m)-TimeOrderShuff(n);
                end
            end
        end
    end
elseif strcmp(input,'delays')
    [NBs NN junk]=size(M); %[NB number X Neuron number X Neuron number
    A=M;
    M=[];
else error('!!!input variable not given correctly'); 
end
%Conversion to Binary delay matrix (unless exists 'NonBinaryWeights')
if strcmp(weights,'Binary')
    A(A>Threshold)=1; 
    A(A<-Threshold)=-1;
    A(A<=Threshold & A>=-Threshold)=0;
    %Matices should be anti-symetric
    AS(AS>Threshold)=1;
    AS(AS<-Threshold)=-1;
    AS(AS<=Threshold & AS>=-Threshold)=0;
    
elseif strcmp(weights,'NonBinary')
else error('!!!weights variable not given correctly'); 
end

hWB=waitbar(0,'Building similarity matrices...');
if averageOverNeurons
    % intialize variables
    R=zeros(NBs,NBs);
    RS=zeros(NBs,NBs);
    RSstd=zeros(NBs,NBs);
    %Build similarity matrix
    for i=1:NBs
        for j=1:NBs
            tmpA=A(i,:,:).*A(j,:,:);
            tmpAS=AS(i,:,:,:).*AS(j,:,:,:);
            
            if ~usePenalty
                tmpA(tmpA<0)=0;
                tmpAS(tmpAS<0)=0;
            end
            %this section should be tested first!!!
                R(i,j)=nansum(tmpA(:))/sum(tmpA(:)~=0 & ~isnan(tmpA(:)));
                RS(i,j)=nansum(tmpAS(:))/sum(tmpAS(:)~=0 & ~isnan(tmpAS(:)));
                RSstd(i,j)=nanstd(nansum(sum(tmpAS,2),3)./sum(nansum(tmpAS~=0  & ~isnan(tmpAS),2),3));
        end
        waitbar(i/NBs,hWB);
    end
else
    % intialize variables
    R=zeros(NBs,NBs,NN);
    RS=zeros(NBs,NBs,NN,nShuffles);
    
    %determine average time order baseline and sort correlation matrices accordingly
    [~,ord]=sort(nanmean(TimeOrder));
    A=A(:,ord,ord);
    AS=AS(:,ord,ord,:);
    
    Mat=~(tril(ones([NN NN]),nearestNeigh) & triu(ones([NN NN]),-nearestNeigh)) & ~eye([NN NN]);
    A(permute(repmat(Mat,[1 1 NBs]),[3 1 2]))=0;
    AS(permute(repmat(Mat,[1 1 NBs nShuffles]),[3 1 2 4]))=0;

    for i=1:NBs
        for j=1:NBs
            tmpA=A(i,:,:).*A(j,:,:);
            tmpAS=AS(i,:,:,:).*AS(j,:,:,:);
            
            if ~usePenalty
                tmpA(tmpA<0)=0;
                tmpAS(tmpAS<0)=0;
            end
            
            R(i,j,:)=nansum(tmpA,2)./sum(tmpA~=0 & ~isnan(tmpA),2);
            RS(i,j,:,:)=nansum(tmpAS,2)./sum(tmpAS~=0 & ~isnan(tmpAS),2);
            %{
            figure;
            h1=subplot(2,4,1);imagesc(conv2(squeeze(M(i,ord,:)),ones(1,5)));
            h2=subplot(2,4,2);imagesc(conv2(squeeze(M(j,ord,:)),ones(1,5)));
            subplot(2,4,3);pcolor(squeeze(A(i,:,:)));
            subplot(2,4,4);pcolor(squeeze(A(j,:,:)));
            subplot(2,4,5);pcolor(squeeze(AS(i,:,:,1)));
            subplot(2,4,6);pcolor(squeeze(AS(j,:,:,1)));
            subplot(2,4,7);pcolor(squeeze(tmpA));
            subplot(2,4,8);pcolor(squeeze(tmpAS(1,:,:,1)));
            colormap(h1,flipud(gray(2)));colormap(h2,flipud(gray(2)));
            
            figure;plot(squeeze(R(i,j,:)),'o-k');hold on;plot(squeeze(RS(i,j,:)),'o-r');plot(squeeze(RSstd(i,j,:)),'o-','color',[0.8 0.5 0.5]);
            %}
        end
        waitbar(i/NBs,hWB);
    end
end
%{
figure;
for i=1:12
    subplot(3,4,i);imagesc(squeeze(M(i,ord,:)));colormap(flipud(gray(2)));
end
%}
close(hWB);