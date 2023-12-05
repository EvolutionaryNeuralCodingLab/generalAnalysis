files=dir('\\storage.laur.corp.brain.mpg.de\Data\Shein-IdelsonMark\Experiments2014\U4_071014\MCRackData\Images3*.mcd');
files={files.name};
dataObj=MCRackRecording(fullfile('\\storage.laur.corp.brain.mpg.de\Data\Shein-IdelsonMark\Experiments2014\U4_071014\MCRackData\',files));

GS=gridSorter(dataObj);
GS.detectionMinimumDetectionInterval=1;
GS=GS.runSorting;