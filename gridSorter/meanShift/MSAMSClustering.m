function clusters=MSAMSClustering(in)
%clustering using automatic meanshift with chooseLandmarks to select the
%initial points
%synthax : clusters=MSAMSClustering(in)
%input : 
%       -in : matrix to clusterize
%output : 
%       - clusters containing index of data points.
[nfeature,nspike]=size(in);    %clustering now
lminfo =choose_landmarks(in,round(nspike/10)+1,struct('frac_random',1,'center',1));
landmarkClust = msams(in,lminfo);
landmarkClust = landmarkClust(lminfo.landmarkAssignment);
[clusters,n_per_group] = agglabel(landmarkClust);
end