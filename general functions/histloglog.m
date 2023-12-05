function [Pdf,T]=histloglog(Isi,res_factor)
% [Pdf,T]=histloglog(Isi);
% calculates the probability density function of the isi as function of t
% on loglog plot. The command loglog(T,Pdf,'.')) will give the desired result.
%res_factor=1.5 is recommended.
%Isi is recommended in msec.
if nargin<2
    res_factor=1.5;
end;
[a b]=size(Isi);
if (a==1) 
    Isi=Isi';
end

%This function looks at positive increments only, to look at both change here
Isi=Isi(find(Isi>=0));
N=length(Isi);

i=3;x(1)=0;x(2)=1;DeltaX(1)=1;DeltaX(2)=1;P(1)=0.5;P(2)=0.5;
while(x(i-1)<1000000)
   x(i)=x(i-1)+DeltaX(i-1);
   DeltaX(i)=round(max(1,10^(res_factor*log10(x(i-1))-2)));
   P(i)=x(i)+DeltaX(i)/2;
   i=i+1;
end

% change 1.1 for 1.5

T=P;
n=histc(Isi,x)';
%[n x]=hist(Isi,x);
%T=x;
Pdf=n./DeltaX/N;
T=T(1:end-1);
Pdf=Pdf(1:end-1);

T(find(Pdf==0))=[];
Pdf(find(Pdf==0))=[];