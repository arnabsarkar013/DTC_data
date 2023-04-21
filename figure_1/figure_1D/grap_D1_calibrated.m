clear all
m=1;

ki=10;kf=60;
 a1=[];
 for k=ki:kf
s1=strcat('grap_',num2str(k),'.csv'); %file path
data1=csvread(s1,42,0);
freq=data1(:,1)./1e6;
Amp1=data1(:,2);
a1=[a1,Amp1];
 end

 a1=10.^((a1-10.*log10(20))/10)./300; 
 a1p=sqrt(a1./(3.8824*10^(11)));


 gate=[0:0.2:10]

figure  
imagesc(gate,freq,a1p);
set(gca,'ColorScale','log')
colormap 'hot'

  set(gca,'YDir','normal');
  set(gca,'XDir','normal');
  xlabel('Pump Voltage(V)')
  ylabel('Frequency(MHz)')