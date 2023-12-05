function [tN,icN]=ManualSort2SortChannel(channels);

eval(['load sort_data_',num2str(channels(1))]);
tN=t;icN=indexchannel;
for i=channels(2:end)
    fprintf('Processing Channel %d\n',i);
    eval(['load sort_data_',num2str(i)]);
    [tN,icN]=MergeSortChannel(tN,icN,t,indexchannel);
end
    