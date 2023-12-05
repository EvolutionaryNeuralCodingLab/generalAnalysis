function groups=groupPairs(E1,E2)
%groupPairs Group sets of link pairs
%   groups=groupPairs(E1,E2) If E1 [N x 1] and E2 [N x 1] are pairs of linked elements (identified by integers), the function returns 
%   cliques merging all connected elements.
%
%	Example:
%	   A=[1 2;2 3;4 1;5 6;7 1];
%	   groups=groupPairs(A(:,1),A(:,2))
%	    
%	See also

% Mark Shein, 26.07.14

nPairs=numel(E1);

groups={};
nGroups=0;
for i=1:nPairs
    %find E1 element in groups or create new group
    currentGroups1=[];
    currentPos1=[];
    for j=1:nGroups
        currentPos1=find(E1(i)==groups{j});
        currentGroups1=j;
        if ~isempty(currentPos1)
            break;
        end
    end
    if isempty(currentPos1)
        nGroups=nGroups+1;
        groups{nGroups}=E1(i);
        currentPos1=1;
        currentGroups1=nGroups;
    end
    
    %find E2 element in groups or create new group
    currentGroups2=[];
    currentPos2=[];
    for j=1:nGroups
        currentPos2=find(E2(i)==groups{j});
        currentGroups2=j;
        if ~isempty(currentPos2)
            break;
        end
    end
    if isempty(currentPos2)
        nGroups=nGroups+1;
        groups{nGroups}=E2(i);
        currentPos2=1;
        currentGroups2=nGroups;
    end
    
    %merge groups
    if currentGroups1~=currentGroups2
        groups{currentGroups1}=[groups{currentGroups1};groups{currentGroups2}];
        groups(currentGroups2)=[];
        nGroups=nGroups-1;
    end
end
