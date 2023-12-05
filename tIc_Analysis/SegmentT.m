function [t_s]=SegmentT(t,hours);
%SegmentT(t,hours) - t in 1/12 [ms];
fprintf('Mininmum value of t is %d hours',min(t)/12000/3600);
segments=ceil(max(t)/12000/3600/hours);
for i=1:segments
    len(i)=length(find(t>(i-1)*hours*3600*12000 & t<i*hours*3600*12000));
end
t_s=zeros(segments,max(len));
for i=1:segments
     tmp_t=t(find(t>(i-1)*hours*3600*12000 & t<i*hours*3600*12000));
     t_s(i,1:length(tmp_t))=tmp_t;
 end;