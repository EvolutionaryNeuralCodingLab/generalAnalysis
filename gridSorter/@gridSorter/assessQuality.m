function [obj, er]=assessQuality(obj)

%determine quantization, repeating in case for backwards compatibility,
%load up spike times and cluster labels, initialize error estimate
%structure
obj.detectionInt2uV=obj.detectionMaxSpikeAmp/2^(obj.detectionNQuantizationBits-1);

load(obj.sortingFileNames.fittingFile,'t','ic','allWaveforms','isNoiseAll','tAll');
load(obj.sortingFileNames.spikeDetectionFile{1},'preSpikeSamplesIntrp');
waveformPeak=preSpikeSamplesIntrp+1;

[nSamples,nNeurons,obj.nCh]=size(allWaveforms);

if ~isempty(ic)
    neuronNames=ic(1:2,:);
else
    neuronNames=[];
end

if isempty(obj.filterObj)
    obj=getHighpassFilter(obj);
end

if ~obj.sortingFileNames.assessQualityExist
    matObj = matfile([obj.sortingDir filesep 'errorEstimates.mat'],'Writable',true);
    matObj.er=struct([]);
    er=struct([]);
else
    load(obj.sortingFileNames.assessQualityFile,'er');
    matObj = matfile([obj.sortingDir filesep 'errorEstimates.mat'],'Writable',true);
    matObj.er=er;
end

fprintf('assessing quality for unit (total %d): ',nNeurons);

for neuron=(size(er,2)+1):nNeurons
    fprintf('%d ',neuron);
    %load up all the surrounding channels,
    %oad all units on surrounding channels
    
    allSpikeFeatures=[];
    allSpikeShapes=[];
    clustLabel=[];
    Th=[];
    
    %for the channels in the small grid surrounding the current neuron's
    %electrode, load up all the spikes into a spike x data point matrix,
    %extract features to approximate nonlinear transofrm of data during
    %clustering
    for c=obj.chPar.surChExtVec{ic(1,neuron)}(obj.chPar.pSurCh{ic(1,neuron)})
        [commonCh,pComN1,pComN2]=intersect(obj.chPar.surChExtVec{ic(1,neuron)}(obj.chPar.pSurCh{ic(1,neuron)}),obj.chPar.surChExtVec{c});
        detectFile=matfile(obj.sortingFileNames.spikeDetectionFile{c});
        skipSpikes=floor(numel(detectFile.spikeTimes)/obj.assessMaxSpikeNumber); %take a subset of spikes to do quality assessment
        currSpikes=detectFile.spikeShapes(:,1:skipSpikes:end,min(pComN2):max(pComN2));
        
        switch obj.featuresFeatureExtractionMethod
            case 'wavelet'
                spikeShapes=double(currSpikes(:,1,pComN2-min(pComN2)+1)) .* obj.detectionInt2uV;
                spikeFeatures=wavedec(spikeShapes,obj.featuresWTdecompositionLevel,obj.featuresSelectedWavelet);
                nCoeffs=numel(spikeFeatures);
                spikeShapes=double(currSpikes(:,:,pComN2-min(pComN2)+1)) .* obj.detectionInt2uV;
                
                spikeFeatures=zeros(size(spikeShapes,2),nCoeffs);
                for j=1:size(spikeShapes,2)
                    spikeFeatures(j,:)=wavedec(spikeShapes(:,j,:),obj.featuresWTdecompositionLevel,obj.featuresSelectedWavelet); %'haar','coif1'
                end
                
                spikeFeatures2=zeros(size(spikeShapes,2),nCoeffs);
                for j=1:size(spikeShapes,2)
                    spikeFeatures2(j,:)=wavedec(permute(spikeShapes(:,j,:),[3 2 1]),obj.featuresWTdecompositionLevel,obj.featuresSelectedWavelet); %'haar','coif1'
                end
                spikeFeatures=[spikeFeatures spikeFeatures2];
                
                
                for j=1:(nCoeffs*2)                               % KS test for coefficient selection
                    thr_dist = std(spikeFeatures(:,j)) * 3;
                    thr_dist_min = mean(spikeFeatures(:,j)) - thr_dist;
                    thr_dist_max = mean(spikeFeatures(:,j)) + thr_dist;
                    aux = spikeFeatures(spikeFeatures(:,j)>thr_dist_min & spikeFeatures(:,j)<thr_dist_max,j);
                    
                    if length(aux) > 10;
                        [ksstat]=test_ks(aux);
                        sd(j)=ksstat;
                    else
                        sd(j)=0;
                    end
                end
                [~,tmp1]=sort(sd(1:nCoeffs),'descend');
                [~,tmp2]=sort(sd(nCoeffs+1:end),'descend');
                spikeFeatures=spikeFeatures(:,[tmp1(1:obj.featuresNWaveletCoeff/2) nCoeffs+tmp2(1:obj.featuresNWaveletCoeff/2)]);
                
                if obj.featuresReduceDimensionsWithPCA
                    [PCAsimMat,spikeFeatures] = princomp(spikeFeatures); %run PCA for visualization purposes
                    spikeFeatures=spikeFeatures(:,1:obj.featuresDimensionReductionPCA);
                end
            case 'PCA' %this option was tested and gives worse results than wavelets
                spikeShapes=double(spikeShapes(:,:,pComN2-min(pComN2)+1)) .* detectionInt2uV;
                [~,spikeFeatures] = princomp(reshape(permute(spikeShapes,[1 3 2]),[nSamples*numel(pComN2-min(pComN2)+1) nSpikes]));
                spikeFeatures=spikeFeatures(1:obj.featuresDimensionReductionPCA,:)';
        end
        
        [time,spks,chan]=size(spikeShapes);
        spikeShapes=reshape(permute(spikeShapes,[1 3 2]),chan*time,spks)';
        allSpikeFeatures=[allSpikeFeatures; spikeFeatures];
        allSpikeShapes=[allSpikeShapes; spikeShapes];
        clustLabel=[clustLabel; tAll{c}(1:skipSpikes:end,1:2)];
        
    end
    
    %find the spikes of the current neuron
    
    iSpikeIndices=find(sum(clustLabel==repmat(ic(1:2,neuron),1,size(clustLabel,1))',2)==2);
    iSpikeFeatures=allSpikeFeatures(iSpikeIndices,:);
    iSpikeShapes=allSpikeShapes(iSpikeIndices,:);
   
 %take error estimates if there are some spikes to work with
if numel(iSpikeIndices)>obj.assessMinSpikeNum
    
    %load up thresholds, top point of peak channel to estimate the percent error estimate due to thresholding
    %detectFile=matfile(obj.sortingFileNames.spikeDetectionFile{ic(1,neuron)});
    Th=detectFile.Th;
    [x loc]=max(abs(mean(iSpikeShapes)));
    spikeAmps=abs(iSpikeShapes(:,loc));
    
    [p,mu,stdev,n,x] = undetected(spikeAmps,mean(Th));
    
    fres3=p; % percent error estimate due to thresholding
    
    %look at overlap with spikes of surrounding units
    confusion=[];
    %allCSpikeFeatures=[];
    allCSpikeLabels=[];
    allCSpikeIndices=[];
    surroundingTally=1;
    channelUnits=unique(clustLabel(:,1));
    for c=1:numel(channelUnits)
        if any(channelUnits(c)==obj.chPar.surChExtVec{ic(1,neuron)}(obj.chPar.pSurCh{ic(1,neuron)})) %don't look at spikes from units not in the surrounding grid
            currChanUnit=clustLabel(find(clustLabel(:,1)==channelUnits(c)),:);
           
            for unit= unique(currChanUnit(:,2))'
                if ~all([currChanUnit(1,1); unit]==ic(1:2,neuron)) %don't double count the unit
                    cSpikeIndices=find(sum(clustLabel==repmat([currChanUnit(1,1); unit],1,size(clustLabel,1))',2)==2);
                    cSpikeFeatures=allSpikeFeatures(cSpikeIndices,:);
                    
                    try
                        confusion(:,:,surroundingTally)=gaussOverlap(iSpikeFeatures,cSpikeFeatures); %fit gaussians to look at overlapping clusters
                    catch
                        confusion(:,:,surroundingTally)=[0 0 ; 0 0];
                    end
                    allCSpikeLabels=[allCSpikeLabels; surroundingTally*ones(numel(cSpikeIndices),1)];
                    allCSpikeIndices=[allCSpikeIndices; cSpikeIndices];
                    surroundingTally=surroundingTally+1;
                end
            end
        end
    end
    
    if ~isempty(confusion)
        fres2(1)=min([max(squeeze(confusion(1,1,:))),1]); fres2(2)=min([max(squeeze(confusion(2,2,:))),1]);  % percent error estimate due to overlapping clusters
    else
        fres2=[0 0];
    end
    er(neuron).fpos=fres2(1)*100; %total false positive rate, in percentage. The sum is super conservative, so I took the max
    er(neuron).fneg=(fres2(2)+fres3)*100; %total false negative rate, in percentage
    % er.fres1=fres1*100; %contamination false positives due to refractory violations
    er(neuron).fres2=fres2*100; %false pos and false neg due to proximity of neighboring clusters
    er(neuron).fres3=fres3*100; %false negatives due to spikes falling below threshold
    er(neuron).performed=1;

else
  er(neuron).fpos=100; %total false positive rate, in percentage. The sum is super conservative, so I took the max
    er(neuron).fneg=100; %total false negative rate, in percentage
    % er.fres1=fres1*100; %contamination false positives due to refractory violations
    er(neuron).fres2=100; %false pos and false neg due to proximity of neighboring clusters
    er(neuron).fres3=100; %false negatives due to spikes falling below threshold
    er(neuron).performed=0;
end


    matObj.er= er;
    
    %% plot stuff for assessment
    if ~exist(['neur_' int2str(neuron) 'qualityAssessment.jpg'])&&numel(iSpikeShapes)>obj.assessMinSpikeNum
        numChan=numel(obj.chPar.pSurCh{ic(1,neuron)});
        numDataPts=size(allWaveforms,1);
        rowSize=surroundingTally+2;
        columnSize=2;
        
        f=figure;
        positionVec=[.05, 1-(2/rowSize), 1/columnSize-0.05, 1/rowSize];
        h=subplot('position',positionVec);
        hold on
        plot(mean(iSpikeShapes))
        
        ylim([min(mean(iSpikeShapes)) max(mean(iSpikeShapes))])
        xlim([0 size(iSpikeShapes,2)])
        for ch=1:numChan
            line([numDataPts*ch numDataPts*ch],[min(mean(iSpikeShapes)) max(mean(iSpikeShapes))],'color','k')
        end
        line([0  numDataPts*numChan],[mean(Th) mean(Th)],'color','r')
        
        set(h,'xtick',[])
        set(h,'xticklabel',[])
        ylabel('uV')
        
        %display the quality assessment and isi plots in top right
        positionVec=[1/columnSize, 1-(2/rowSize), 1/(columnSize*2)-0.05, 1/rowSize];
        h=subplot('position',positionVec);
        hold on
        text(.010,.8,['false pos: ' int2str(er(neuron).fpos) '%'],'Units','normalized')
        text(.010,.5,['false neg: ' int2str(er(neuron).fneg) '%'],'Units','normalized')
        text(.010,.2,['thresh errors: ' int2str(er(neuron).fres3) '%'],'Units','normalized')
        % text(.8,.8,[int2str(size(iSpikeShapes,1)) ' spikes'],'Units','normalized')%
        set(h,'xtick',[])
        set(h,'xticklabel',[])
        set(h,'ytick',[])
        set(h,'yticklabel',[])
        
        positionVec=[1.5*(1/columnSize)+0.05, 1-(2/rowSize), 1/(columnSize*2)-0.1, 1/rowSize];
        h=subplot('position',positionVec);
        hold on
        plot(0:0.01:1,histc(diff(t(ic(3,neuron):ic(4,neuron))),0:0.01:1))
        title('ISI histogram (0:1 ms)')
        ylabel('spikes')
        set(h,'xtick',[])
        set(h,'xticklabel',[])
        
        if surroundingTally>=1
            for unit=1:(surroundingTally-1)
                cSpikeIndices=allCSpikeIndices(find(allCSpikeLabels==unit));
                cSpikeShapes=allSpikeShapes(cSpikeIndices,:);
                
                positionVec=[.05, 1-((2+unit)/rowSize), 1/columnSize-0.05, 1/rowSize];
                h=subplot('position',positionVec);
                plot(mean(cSpikeShapes),'r')
                for c=1:numChan
                    line([numDataPts*c numDataPts*c],[min(mean(iSpikeShapes)) max(mean(iSpikeShapes))],'color','k')
                end
                ylim([min(mean(iSpikeShapes)) max(mean(iSpikeShapes))])
                xlim([0 numDataPts*numChan])
                set(h,'xtick',[])
                set(h,'xticklabel',[])
                set(h,'ytick',[])
                set(h,'yticklabel',[])
                
                positionVec=[1/columnSize, 1-((2+unit)/rowSize), 1/columnSize-0.05, 1/rowSize];
                h=subplot('position',positionVec);
                hold on
                label=[ones(1,size(cSpikeShapes,1))+1 ones(1,size(iSpikeShapes,1))];
                relSpikes=[allSpikeFeatures(cSpikeIndices,:); iSpikeFeatures];
                
                try
                    [redDimSpikes, x1] = lda(relSpikes, label, 2);
                    scatter( redDimSpikes(:,1),redDimSpikes(:,2),20,label,'filled')
                    
                catch
                end
                y1=ylim;
                x1=xlim;
                line([x1(1) x1(2)],[y1(2) y1(2)],'color','k')
                % text(.8,.8,[int2str(size(cSpikeShapes,1)) ' spikes'],'Units','normalized')
                set(h,'xtick',[])
                set(h,'xticklabel',[])
                set(h,'ytick',[])
                set(h,'yticklabel',[])
            end
        end
        set(f,'units','normalized','outerposition',[.1 .1 .8 .8])
        printFile=[obj.sortingDir filesep 'neur_' int2str(neuron) 'qualityAssessment'];
        set(f,'PaperPositionMode','auto');
        print(printFile,'-djpeg','-r300');
        close all
    end
    
end

for x=1:size(er,2)
    fpos(x)=er(x).fpos;
    fneg(x)=er(x).fneg;
    thresh(x)=er(x).fres3;
end

figure
hist(thresh,30)
title('errors due to thresholding')
xlabel('error rate (%)')
ylabel('count')
xlim([0 100])
figure
hist(fpos,30)
title('errors due to cluster overlap false positives')
xlabel('error rate (%)')
ylabel('count')
xlim([0 100])
figure
hist(fneg-thresh,30)
title('errors due to cluster overlap false negatives')
xlabel('error rate (%)')
ylabel('count')
xlim([0 100])



%

    function [p,mu,stdev,n,x] = undetected(w,threshes)
        % UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/201
        %modified by Sam Reiter 7.4.2015
        
        % Output:
        %  p            - estimate of probability that a spike is missing because it didn't reach threshhold
        %  mu           - mean estimated for gaussian fit
        %  stdev        - standard deviation estimated for gaussian fit
        %  n            - bin counts for histogram used to fit Gaussian
        %  x            - bin centers for histogram used to fit Gaussian
        %
        bins=75;
        % normalize all waveforms by threshold
        w = abs(w);
        w=w./ abs(threshes);
        
        
        % create the histogram values
        global_max = max(w);
        mylims = linspace( 1,global_max,bins+1);
        x = mylims +  (mylims(2) - mylims(1))/2;
        n = histc( w,mylims );
        
        % fit the histogram with a cutoff gaussian
        m = mode_guesser(w, .05);    % use mode instead of mean, since tail might be cut off
        [stdev,mu] = stdev_guesser(w, n, x, m); % fit the standard deviation as well
        
        % Now make an estimate of how many spikes are missing, given the Gaussian and the cutoff
        p = normcdf( 1,mu,stdev);
        
        % attempt to keep values negative if all threshold values were negative
        if all( threshes < 0 )
            mu = -mu;
            x = -x;
        end
        
    end



% fit the standard deviation to the histogram by looking for an accurate
% match over a range of possible values

    function [stdev,m] = stdev_guesser( thresh_val, n, x, m)
        % initial guess is juts the RMS of just the values below the mean
        init = sqrt( mean( (m-thresh_val(thresh_val>=m)).^2  ) );
        
        % try 20 values, within a factor of 2 of the initial guess
        num = 20;
        st_guesses = linspace( init/2, init*2, num );
        m_guesses  = linspace( m-init,max(m+init,1),num);
        for j = 1:length(m_guesses)
            for k = 1:length(st_guesses)
                b = normpdf(x,m_guesses(j),st_guesses(k));
                b = b *sum(n) / sum(b);
                error(j,k) = sum(abs(b(:)-n(:)));
            end
        end
        
        % which one has the least error?
        [val,pos] = min(error(:));
        jpos = mod( pos, num ); if jpos == 0, jpos = num; end
        kpos = ceil(pos/num);
        stdev = st_guesses(kpos);
        
        % refine mode estimate
        m     = m_guesses(jpos);
        
    end


    function confusion = gaussOverlap( w1, w2 )
        % UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
        %modified by Sam Reiter 7.4.2015
        
        % Input:
        %   waveforms1  - [Event x Sample ] waveforms of 1st cluster
        %   waveforms2  - [Event x Sample ] waveforms of 2nd cluster
        %      If your waveform data is in the form [Event X Sample X Channels ],
        %      you should call this function as
        %            confusion = gaussian_overlap( w1(:,:), w2(:,:) );
        %      This will concatenate the waveforms from different channels.
        %
        % Output:
        %   C   - a confusion matrix
        %   C(1,1) - False positive fraction in cluster 1 (waveforms of neuron 2 that were assigned to neuron 1)
        %   C(1,2) - False negative fraction in cluster 1 (waveforms of neuron 1 that were assigned to neuron 2)
        %   C(2,1) - False negative fraction in cluster 2
        %   C(2,2) - False positive fraction in cluster 2
        %
        
        % fit 2 multivariate gaussians, use observed parameters to initialize
        params.mu = [ mean(w1(:,:),1); mean(w2(:,:),1)];
        params.Sigma(:,:,1) = cov(w1(:,:));
        params.Sigma(:,:,2) = cov(w2(:,:));
        
        % disp('Fitting 2-Gaussian Mixture Model ...')
        gmfit = gmdistribution.fit([w1(:,:); w2(:,:)],2,'Start',params);
        
        % get posteriors
        %disp('Calculating Confusion matrix ...')
        pr1 = gmfit.posterior(w1);
        pr2 = gmfit.posterior(w2);
        
        % in the unlikely case that the cluster identities were flipped during the fitting procedure, flip them back
        if mean(pr1(:,1)) + mean(pr2(:,2)) < 1
            pr1 = pr1(:,[2 1]);
            pr2 = pr2(:,[2 1]);
        end
        
        % create confusion matrix
        confusion(1,1) = mean(pr1(:,2));    % probability that a member of 1 is false
        confusion(1,2) = sum(pr2(:,1))/size(w1,1); % relative proportion of spikes that were placed in cluster 2 by mistake
        confusion(2,1) = sum(pr1(:,2))/size(w2,1); %  relative proportion of spikes that were placed in cluster 1 by mistake
        confusion(2,2) = mean(pr2(:,1));   % probability that a member of 2 was really from 1
    end



    function m = mode_guesser(x,p)
        
        % UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/09/2010
        %
        % mode_guesser - guess mode of the
        %
        % Usage:
        %    m = mode_guesser(x,p)
        %
        % Description:
        %   Guesses mode by looking for location where the data is most tightly
        % distributed.  This is accomplished by sorting the vector x and
        % looking for the p*100 percentile range of data with the least range.
        %
        % Input:
        %   x - [1 x M] vector of scalars
        %
        % Option input:
        %   p - proportion of data to use in guessing the mode, defaults to 0.1
        %
        % Output:
        %   m - guessed value of mode
        %
        
        %check for whether p is specified
        if nargin < 2, p = .1; end
        
        % determine how many samples is p proportion of the data
        num_samples = length(x);
        shift = round( num_samples * p );
        
        % find the range of the most tightly distributed data
        x = sort(x);
        [val,m_spot] = min( x(shift+1:end) - x(1:end-shift) );
        
        % use median of the tightest range as the guess of the mode
        m = x( round(m_spot + (shift/2)) );
    end



    function [mappedX, mapping] = lda(X, labels, no_dims)
        %LDA Perform the LDA algorithm
        %
        %   [mappedX, mapping] = lda(X, labels, no_dims)
        %
        % The function runs LDA on a set of datapoints X. The variable
        % no_dims sets the number of dimensions of the feature points in the
        % embedded feature space (no_dims >= 1, default = 2). The maximum number
        % for no_dims is the number of classes in your data minus 1.
        % The function returns the coordinates of the low-dimensional data in
        % mappedX. Furthermore, it returns information on the mapping in mapping.
        %
        %
        
        % This file is part of the Matlab Toolbox for Dimensionality Reduction.
        % The toolbox can be obtained from http://homepage.tudelft.nl/19j49
        % You are free to use, change, or redistribute this code in any way you
        % want for non-commercial purposes. However, it is appreciated if you
        % maintain the name of the original author.
        %
        % (C) Laurens van der Maaten, Delft University of Technology
        
        
        if ~exist('no_dims', 'var') || isempty(no_dims)
            no_dims = 2;
        end
        
        % Make sure data is zero mean
        mapping.mean = mean(X, 1);
        X = bsxfun(@minus, X, mapping.mean);
        
        % Make sure labels are nice
        [classes, bar, labels] = unique(labels);
        nc = length(classes);
        
        % Intialize Sw
        Sw = zeros(size(X, 2), size(X, 2));
        
        % Compute total covariance matrix
        St = cov(X);
        
        % Sum over classes
        for i=1:nc
            
            % Get all instances with class i
            cur_X = X(labels == i,:);
            
            % Update within-class scatter
            C = cov(cur_X);
            p = size(cur_X, 1) / (length(labels) - 1);
            Sw = Sw + (p * C);
        end
        
        % Compute between class scatter
        Sb       = St - Sw;
        Sb(isnan(Sb)) = 0; Sw(isnan(Sw)) = 0;
        Sb(isinf(Sb)) = 0; Sw(isinf(Sw)) = 0;
        
        % Make sure not to embed in too high dimension
        if nc < no_dims
            no_dims = nc;
            warning(['Target dimensionality reduced to ' num2str(no_dims) '.']);
        end
        
        % Perform eigendecomposition of inv(Sw)*Sb
        [M, lambda] = eig(Sb, Sw);
        
        % Sort eigenvalues and eigenvectors in descending order
        lambda(isnan(lambda)) = 0;
        [lambda, ind] = sort(diag(lambda), 'descend');
        M = M(:,ind(1:min([no_dims size(M, 2)])));
        
        % Compute mapped data
        mappedX = X * M;
        
        % Store mapping for the out-of-sample extension
        mapping.M = M;
        mapping.val = lambda;
        
        
    end
end
