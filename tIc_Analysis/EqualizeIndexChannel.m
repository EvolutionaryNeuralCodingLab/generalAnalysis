function [newT1,newIndexchannel1,newT2,newIndexchannel2]=EqualizeIndexChannel(t1,indexchannel1,t2,indexchannel2);
%function [newT1,newIndexchannel1,newT2,newIndexchannel2]=EqualizeIndexChannel(t1,indexchannel1,t2,indexchannel2);
%a functioin that returns newT1,2 and newIndexchannel1,2 of common neurons only

newT1=t1;
newT2=t2;
newIndexchannel1=indexchannel1;
newIndexchannel2=indexchannel2;
for i=1:length(indexchannel1),
    if isempty(find(indexchannel2(1,:)==indexchannel1(1,i) & indexchannel2(2,:)==indexchannel1(2,i))),
        [newT1,newIndexchannel1]=RemoveNeurons(newT1,newIndexchannel1,indexchannel1(1:2,i));
    end
end
for i=1:length(indexchannel2),
    if isempty(find(indexchannel1(1,:)==indexchannel2(1,i) & indexchannel1(2,:)==indexchannel2(2,i))),
        [newT2,newIndexchannel2]=RemoveNeurons(newT2,newIndexchannel2,indexchannel2(1:2,i));
    end
end
