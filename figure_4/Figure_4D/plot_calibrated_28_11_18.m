clear all
m=1;
 ki=10;kf=135;

 a=[];
 a1=[];
 a2=[];
 for k=ki:kf
 
s=strcat('det_',num2str(k),'.csv'); %file path

data1=csvread(s,42,0);

freq=data1(:,1)./1e6;
Amp=data1(:,2);

Amp1=max(Amp(178:214));
a1=[a1,Amp1];
Amp2=max(Amp(220:250));
a2=[a2,Amp2];
a=[a,Amp];
  
end
a=10.^((a-10.*log10(20))/10)./300;
a1=10.^((a1-10.*log10(20))/10)./300;
a2=10.^((a2-10.*log10(20))/10)./300;
% a=sqrt(a./(8.5984*1e16));
a=sqrt(a./(3.8824*10^(11)));
a1=sqrt(a1./(3.8824*10^(11)));
a2=sqrt(a2./(3.8824*10^(11)));



  gate=[5.3:0.004:5.8]

imagesc(gate,freq,(a));
set(gca,'ColorScale','log')
   colormap hot
  set(gca,'YDir','normal');
  ylabel('Frequency(MHz)')
  xlabel('Detuning(MHz)')