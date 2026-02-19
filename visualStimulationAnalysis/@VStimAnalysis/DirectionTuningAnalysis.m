%%%%%%% tuning analysis

function fig = DirectionTuningAnalysis(expList, Stims2Comp,params)

arguments
    expList  (1,:) double  %%Number of experiment from excel list
    Stims2Comp cell %% Comparison order {'MB','RG','MBR'} would select neurons responsive to moving ball and 
    % compare this neurons responses to other stimuli. 
    params.threshold = 0.05;
    params.diffResp = false;
    params.overwrite = false;
    params.StimsPresent = {'MB','RG'}; %assumes that at least moving ball is present
    params.StimsNotPresent = {};
    params.StimsToCompare = {}; %Select 2 stims to compare scatter plots (default: 1st and 2nd stim are compared from the Stims2Comp cell array)
    params.overwriteResponse = false;
    params.RespDurationWin = 100; %same as default
    params.shuffles = 2000; %same as default
    params.ignoreNonSignif = false; %when comparing first stim, ignore neurons non responsive to other stim
    params.EachStimSignif = false; %resposnive neurons for each stim are selected (default: responsive neurons of first stime are selected)
end


cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

GoodRecordingsPV =[49:54];

recordings = GoodRecordingsPV;

takeMedian = 1;

tuningCurve = cell(1,numel(recordings));
SEM = cell(1,numel(recordings));
DSI = cell(1,numel(recordings));
OSI =cell(1,numel(recordings));
Theta = cell(1,numel(recordings));
preferDir = cell(1,numel(recordings));

MBc =[];
insertionID = [];
i =1;
animal =0;
for ex =  GoodRecordingsPV
    %%%%%%%%%%%% Load data and data paremeters

    %1. Load NP class
    path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    try %%In case it is not run in Vstim computer, which has drives mapped differently
        cd(pathE)
    catch
        try
            originP = cell2mat(extractBetween(path,"\\","\Large_scale"));
            if strcmp(originP,'sil3\data')
                path = replaceBetween(path,"","\Large_scale","W:");
            else
                path = replaceBetween(path,"","\Large_scale","Y:");
            end

            cd(path)
        catch

            if strcmp(originP,'sil3\data')
                path = replaceBetween(path,"","\Large_scale","\\sil3\data");
            else
                path = replaceBetween(path,"","\Large_scale","\\sil1\data");
            end
            cd(path)

        end
    end
    NP = NPAPRecording(path);

    N_bootstrap =1000;
    cd(NP.recordingDir)
    respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
    sign = 0.05; %%%Significance level used to calculate receptive fields
    respU = find(respNeuronsMB<sign);
        
    tempTC = load('tuningCurveAllOffsets').tuningCurve;;

    if size(tempTC,2) > 4
        tempTC = tempTC(:,1:2:size(tempTC,2)); %%Select only up, left, down, and right (0,90,180,270) if the recording has more than 4 directions
    end


    tuningCurve{i} = tempTC;


    DSIt= load(sprintf('Direction-Selectivity-Index-%s',NP.recordingName)).DSI;
    DSIt = DSIt(respNeuronsMB<sign);
    DSI{i} = DSIt;
    
    OSIt = load(sprintf('Orientation-Tuning-Index-%s',NP.recordingName)).L;
    OSIt = OSIt(respNeuronsMB<sign);
    OSI{i} = OSIt;
    
    Thetat = load(sprintf('Angle-prefer-%s',NP.recordingName)).pI;
    Thetat = Thetat(respNeuronsMB<sign);
    Theta{i} = Thetat;
    [maxVal preferD] = max(tempTC(respNeuronsMB<sign,:),[],2);

    preferDir{i} = preferD;

     %%%Add color vectors accordying to animals
    if string(data.Animal_ID{ex}) ~= string(data.Animal_ID{ex-1}) %wont work if you start with the first animal (noisy animal)
        animal = animal+1;
        animalName{animal} = data.Animal_ID{ex};
    end

    MBc = [MBc ...
            zeros(1,length(DSI{i}))+ animal]; %%%Add animal number for colors

    insertionID = [insertionID repmat(i,1,numel(DSIt))];

    i=i+1;

end


%% Colors per insertion
in =i;
% Number of colors
nColors = in-1;

% Generate hues evenly spaced around the color wheel
H = linspace(0, 1, nColors+1); 
H = H(1:end-1); % Remove last to avoid duplication at 1 and 0

% Keep saturation low for pale colors
S = 0.5 * ones(1, nColors); 

% Keep value high for brightness
V = 0.7 * ones(1, nColors); 

% Convert HSV to RGB
colorMatrix = hsv2rgb([H' S' V']);

% Display the colors
figure;
colormap(colorMatrix);
c = colorbar;
caxis([0 1])
c.Ticks = linspace(1/(in-1),1,in-1)-(1/(in-1))/2; 
c.TickLabels = 1:in-1;
c.Limits
title('Pale Differentiated Colors');
legend()

colors = colorMatrix(insertionID,:);


%% plot histogram with a histogram per direction that shows the distribution of OSI and DSI values

% cats = categorical([ones(1,length(cell2mat(DSI'))) ones(1,length(cell2mat(DSI')))+1 ones(1,length(cell2mat(DSI')))+2])';
% 
% fig = figure;
% tiledlayout(2,1, "TileSpacing","loose","Padding","loose")
% 
% %%%%%%%%%%%
% nexttile %%% all
% 
% swarmchart(cats, [cell2mat(DSI');cell2mat(OSI');cell2mat(DSI')-cell2mat(OSI')], 10, [colors;colors;colors], 'filled','MarkerFaceAlpha',0.5)
% 
% set(gcf,'Color','w');%
% xticklabels({'DSI','OSI','DSI-OSI'})
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 7; % Set font size of x-axis tick labels
% ylabel('Tuning strength')
% yline(0,'LineWidth',1,'Color','k')

vDSI = cell2mat(DSI');

vOSI = cell2mat(OSI');
% 

% %%%%%%%%%%
% nexttile %%% dsi > 0.7
% 
% thres1 = 0.7;

% 
% diffDSI_OSI = vDSI(vDSI>thres1)-vOSI(vDSI>thres1);
% 
% val = [vDSI(vDSI>thres1);vOSI(vDSI>thres1);diffDSI_OSI];
% 
% cats = categorical([ones(1,sum(vDSI>thres1)) ones(1,sum(vDSI>thres1))+1 ones(1,sum(vDSI>thres1))+2])';
% 
% colors = [MBc(vDSI>thres1)'; MBc(vDSI>thres1)';MBc(vDSI>thres1)'];
% 
% swarmchart(cats, val, 10,'filled','MarkerFaceAlpha',0.8)
% 
% yline(thres1,'LineWidth',2,'Label',string(0.7),'LabelHorizontalAlignment','right','FontSize',7)
% xticklabels({'DSI > 0.7','OSI','DSI-OSI'})
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 7; % Set font size of x-axis tick labels
% ylabel('Tuning strength')


% %%%%%%%
% nexttile %%% osi > 0.7
% 
% 
% thres1 = 0.7;
% 
% 
% diffOSI_DSI = vOSI(vOSI>thres1)-vDSI(vOSI>thres1);
% 
% val = [vDSI(vOSI>thres1);vOSI(vOSI>thres1);diffOSI_DSI];
% 
% cats = categorical([ones(1,sum(vOSI>thres1)) ones(1,sum(vOSI>thres1))+1 ones(1,sum(vOSI>thres1))+2])';
% 
% colors = [MBc(vOSI>thres1)'; MBc(vOSI>thres1)';MBc(vOSI>thres1)'];
% 
% swarmchart(cats, val, 10,'filled','MarkerFaceAlpha',0.8)
% 
% yline(thres1,'LineWidth',2,'Label',string(0.7),'LabelHorizontalAlignment','right','FontSize',7)

% set(gcf,'Color','w');%
% % xticklabels({'DSI','OSI > 0.7','OSI-DSI'})
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 7; % Set font size of x-axis tick labels



ylabel('Tuning strength')
ylim([-0.1 1])
cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
fig = figure;
scatter(vDSI,vOSI,7, colors,"filled",'MarkerFaceAlpha',0.7);
xlabel('DSI')
ylabel('OSI')
axis equal
xlim([0 1])
ylim([0 1])
hold on
plot([0,1],[0,1],'LineWidth',1,'Color','k')
fig.Position = [1141         261         214         140];

print(fig, 'tuningIndexesMovBall.pdf', '-dpdf', '-r300', '-vector');

%% Histogram of prefered angles

%%Create a round histogram divided into 4 direction. Within each direction,
%%DSI values are sorted

angles = cell2mat(preferDir');
DSIv = cell2mat(DSI');
OSIv = cell2mat(OSI');

figure;swarmchart(categorical(angles),OSIv,10, colors,'filled','MarkerFaceAlpha',0.8);

set(gcf,'Color','w');%
ax = gca; % Get current axis
ax.XAxis.FontSize = 7; % Set font size of x-axis tick labels
ylabel('DSI')
yline(0,'LineWidth',1,'Color','k')

fig.Position = [1269         362         189         305];

print(gcf, 'tuningIndexesMovBallPerAngle.pdf', '-dpdf', '-r300', '-vector');


end








