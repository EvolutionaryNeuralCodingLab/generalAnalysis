function [gc,f]=projectionMerge(spikeFeatures,initIdx,varargin)
%merge clusters in a group that have to be merge based on the residuals between spikes and templates.
%synthax : [gc,f]=projectionMerge(spikeFeatures,initIdx,varargin)
%input:
%           - spikeFeatures : spike features ()
%           - initIdx : the clusters index of every spike
%               vararin - 'property','value'
%output :
%           - gc: a binary matrix with ones indicating a necessary merge
%           - f: a figure handle for the generated plot

%default variables
obj.clusteringMinNSpikesCluster=10;
obj.clusteringSTDMergeFac=2;
obj.clusteringMergeThreshold=0.18;
obj.clusteringPlotProjection=1;

%Collects all options
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

%find robust cluster centers
nClustersIn=numel(unique(initIdx));
for k=1:nClustersIn
    cent(k,:)=median(spikeFeatures(initIdx==k,:));
end

if  nClustersIn<=1
    gc=1;
    f=[];
    return;
end

if obj.clusteringPlotProjection
    f=figure('Position',[50 50 1400 900]);
else
    f=[];
end

D=zeros(nClustersIn);
groups=mat2cell(1:nClustersIn,1,ones(1,nClustersIn));
gc=1:nClustersIn; % a group is assigned to every cluster
for k=2:nClustersIn
    for m=1:k-1
        v=(cent(k,:)-cent(m,:));
        p1=projectionND(v,spikeFeatures(initIdx==k,:));
        p2=projectionND(v,spikeFeatures(initIdx==m,:));
        
        pCent=projectionND(v,[cent(k,:);cent(m,:)]);%for plotting purpuses
        
        if numel(p1)<obj.clusteringMinNSpikesCluster || numel(p2)<obj.clusteringMinNSpikesCluster
            D(k,m)=0;
            std_p1=[];
            std_p2=[];
            v=[];
        else
            
            mp1=median(p1);
            mp2=median(p2);
            
            std_p1=median(abs(  p1-mp1  ),2) / 0.6745;
            std_p2=median(abs(  p2-mp2  ),2) / 0.6745;
            
            %std_p1=std(p1);
            %std_p2=std(p1);
            nV=numel(p1)+numel(p2);
            
            s=sign(mp1-mp2);
            %edges=[(mp2-std_p2*s):(s*10/nV*abs(mp1-mp2)):(mp1+std_p1*s)];
            edges=(mp2-std_p2*s):(((mp1+std_p1*s)-(mp2-std_p2*s))/(round(log(nV))/10)/20):(mp1+std_p1*s);
            %eges must be divided in a way that preserve the extreme edges on both sides
            n1=histc(p1,edges);
            n2=histc(p2,edges);
            n=n2(1:end-1)-n1(1:end-1); %eliminate edges that sum over
            %edges=edges(1:end-1);
            
            firstCross=find(n(2:end)<=0 & n(1:end-1)>0,1,'first')+1;
            secondCross=find(n(1:end-1)>=0 & n(2:end)<0,1,'last');
            intersection=(edges(firstCross) + edges(secondCross))/2;
            %D(k,m)=max(sum(p2>(intersection-std_p2/obj.clusteringSTDMergeFac))/sum(n2),sum(p1<(intersection+std_p1/obj.clusteringSTDMergeFac))/sum(n1));
            %D(k,m)=sum([p1 p2]<(intersection+std_p2/obj.clusteringSTDMergeFac) & [p1 p2]>(intersection-std_p1/obj.clusteringSTDMergeFac))/sum([n1 n2]);
            if isempty(intersection) %one cluster is contained within the other
                D(k,m)=1;
            else
                D(k,m)=(  sum(p2>(intersection-std_p2/obj.clusteringSTDMergeFac))  +  sum(p1<(intersection+std_p1/obj.clusteringSTDMergeFac))   ) /sum([n1 n2] );
            end
            %figure;plot([p1 p2]);hold on;plot(ones(1,numel(edges)),edges,'or');line([1 sum([n1 n2])],[intersection intersection],'color','k','LineWidth',2);
        end
        %D(k,m)=(sqrt(sum(v.^2))/sqrt(std_p1^2 + std_p2.^2));
        %D(k,m)=(sqrt(sum(v.^2))/(std_p + std_p2);
        %D(k,m)=(1+skewness([p1 p2])^2)/(kurtosis([p1 p2])+3);
        %D(k,m) = kstest2(p1,p2,[],0.0.05);
        
        if D(k,m)>obj.clusteringMergeThreshold
            %find in which group k is and add all is group to m
            groupOfK=gc(k);
            groupOfM=gc(m);
            gc(gc==groupOfK)=groupOfM;
        end
        
        if obj.clusteringPlotProjection
            
            subaxis(nClustersIn,nClustersIn,(m-1)*nClustersIn+k,'Spacing', 0.001, 'Padding', 0.001, 'Margin', 0.001);
            
            edges=[min([p1 p2]):((max([p1 p2])-min([p1 p2]))/30):max([p1 p2])]; %edges different from before
            n1=hist(p1,edges);
            n2=hist(p2,edges);
            bar(edges,[n1;n2]',1,'stacked');
            axis tight;
            set(gca,'XTickLabel',[],'YTickLabel',[]);
            
            subaxis(nClustersIn,nClustersIn,(k-1)*nClustersIn+m,'Spacing', 0.001, 'Padding', 0.001, 'Margin', 0.001);
            strTxt={['F=' num2str(D(k,m))],['s1=' num2str(std_p1)],['s2=' num2str(std_p2)],['D=' num2str(sqrt(sum(v.^2)))]};
            text(0,0.5,strTxt);
            axis off;
        end
        
    end
end

    function p=projectionND(v,d)
        %Calculate projection between a vector and a set of dots in multi dimensional space
        %v = [1 x N] - vector
        %d = [M X N] - M dot locations
        
        %calculate the cos angle between vector and dots
        cosAng=v*d'./(sqrt(sum(v.^2))*sqrt(sum(d'.^2)));
        p=cosAng.*sqrt(sum(d'.^2));
    end
end