function [tnew,indexchannelnew]=InsertNeuron(t,indexchannel,neuronT,neuronIndex);
%function [tnew,indexchannelnew]=InsertNeuron(t,indexchannel,neuronT,neuronIndex);
%a fucntion that adds a new neuron to t, indexchannel;

indexchannelnew=[indexchannel,[ neuronIndex;length(t)+1;0]];
tnew=[t,neuronT];
indexchannelnew(4,end)=length(tnew);

    