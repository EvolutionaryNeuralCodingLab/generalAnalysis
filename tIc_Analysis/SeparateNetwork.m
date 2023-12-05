function [t_a,ic_a,t_b,ic_b]=separateNetwork(t,ic,varargin)
%[t_a,ic_a,t_b,ic_b]=separateNetwork(t,ic,'500um');
%
% Function purpose :    separates a network into 2 parts A and B
% Function recives :    t - firing timings
%                       ic - indexchannel
%                       varargin:
%                       'separation'='500um','1-50,51-100','manual'
% Function give back :  t_a,t_b - the timings of both separated networks
%                       ic_a,ic_b - the index channel of both separated networks
%
% Last updated : 10/10/14

%default variables
separation='1-50,51-100';
outputMessages=false;

%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
    end
    return;
end

%collect input variables
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

if strcmp(separation,'500um')
    side_a=[1 2 3 7 8 9 10 15 16 17 18 23 24 25 26 31 32 33 34 39 40 41 42 47 48 49 50 55 56 57];
    side_b=[4 5 6 11 12 13 14 19 20 21 22 27 28 29 30 35 36 37 38 43 44 45 46 51 52 53 54 58 59 60];
elseif strcmp(separation,'1-50,51-100')
    side_a=1:50;
    side_b=51:100;
end

count=1;
a=[];b=[];
t_a=t;t_b=t;ic_a=ic;ic_b=ic;
if strcmp(separation,'manual')
    fprintf('For the following neurons write a for network a and b for b \n');
    while size(ic,2) >= count
        fprintf('channel %d neuron %d -',ic(1,count),ic(2,count));
        answ=input(' ','s');
        if answ=='a'
            a=[a ic(1:2,count)];
            count=count+1;
        else
            if answ=='b'
                b=[b ic(1:2,count)];
                count=count+1;
            else
                fprintf('\n Your answer is not a or b,try again ! \n');
            end
        end
    end
else
    for i=1:size(ic,2)
        if ~isempty(find(side_a==ic(1,i)))
            a=[a ic(1:2,i)];
        else
            b=[b ic(1:2,i)];
        end
    end
end
if outputMessages
    fprintf('Removing the Neurons from A:\n');
end
for i=1:size(b,2)
    [t_a,ic_a]=removeNeurons(t_a,ic_a,b(:,i));
end
if outputMessages
    fprintf('Removing the Neurons from B:\n');
end
for i=1:size(a,2)
    [t_b,ic_b]=removeNeurons(t_b,ic_b,a(:,i));
end

