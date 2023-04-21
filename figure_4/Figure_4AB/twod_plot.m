%Frquency mixing in graphene resonator
clear all
m=1;

  ki=10;kf=110;

%  a=[];
 a1=[];
 a2=[];
% a=zeros
 for k=ki:kf

s1=strcat('grap_',num2str(k),'.csv'); %file path
% s2=strcat('DC_scan_D2\grap_',num2str(k),'.csv'); %file path

data1=csvread(s1,42,0);
% data2=csvread(s2,42,0);

freq=data1(:,1)./1e6;

Amp1=data1(:,2);
% Amp2=data2(:,2);

% Amp1=max(Amp(178:214));
% a1=[a1,Amp1];
% Amp2=max(Amp(220:250));
% a2=[a2,Amp2];
a1=[a1,Amp1];
% a2=[a2,Amp2];
 end

 a1=10.^((a1-10.*log10(20))/10)./10; 
%  a2=10.^((a2-10.*log10(20))/10)./300; 
%  a=a1+a2;
%  a=sqrt(a./(3.8824*10^(11)));
 a1p=sqrt(a1./(3.8824*10^(11)));
%  a2p=sqrt(a2./(3.8824*10^(11)));

 gate=[10:-0.1:0]

figure(2)
imagesc(gate,freq,(a1p));
set(gca,'ColorScale','log')
colormap 'hot'

  set(gca,'YDir','normal');
  set(gca,'XDir','normal');
  xlabel('Pump Voltage(V)')
  ylabel('Frequency(MHz)')
% set(gca,'zscale','log')