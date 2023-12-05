% [MC]=ConvBurstMatrix(M,kern,shape);
% Function purpose : convolutes multiple one dimensional vectors
%
% Recives :             M [mxn] - a matrix with n channels, each with m samples
%                                   kern [1xn] - the kernel for convolution 
%                                   shape - returns a subsection of the N-dimensional convolution, as specified by the shape parameter:
%                                       'full' - Returns the full N-dimensional convolution (default).
%                                       'same' - Returns the central part of the result that is the same size as A.
%                                       'valid' - Returns only those parts of the convolution that can be computed without assuming that the array A is zero-padded. 
%                                                                                                                                                        
% Function give back :  MC - convoluted matrix
% Recomended usage  : [MC]=ConvBurstMatrix(M,fspecial('gaussian',[1 10],3),'same');
% Last updated : 28/03/11
function [MC]=ConvBurstMatrix(M,kern,shape)

%Chech input parameter and set default values
if size(kern,1)>1
    error('The input kernel is not one dimensional');
end
if (nargin < 3)
  shape = 'full';
end

%Initialize output variable
[NBursts NNeurons Width]=size(M);
NKern=size(kern,2);
if (strcmp(shape,'same'))
    MC=zeros(NBursts,NNeurons,Width);
elseif (strcmp(shape,'full'))
    MC=zeros(NBursts,NNeurons,Width+NKern-1);
else
    error('The only supported shapes are ''full'' and ''same''');
end

if NBursts==1
    MC(1,:,:)=convn(squeeze(M(1,:,:))',kern,shape);
else
    for i=1:NBursts
        MC(i,:,:)=convn(squeeze(M(i,:,:)),kern,shape);
    end
end

%Calculate convolution
