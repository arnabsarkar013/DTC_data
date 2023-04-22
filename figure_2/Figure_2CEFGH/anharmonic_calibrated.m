clear all
m=1;
  ki=10;kf=110;
 a=[];
 a1=[];
 a2=[];
 for k=ki:kf

s=strcat('anharmonic_D1\grap_',num2str(k),'.csv'); %file path

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
a=sqrt(a./(3.8824*10^(11)));
a1=sqrt(a1./(3.8824*10^(11)));
a2=sqrt(a2./(3.8824*10^(11)));


gate=[10:-0.1:0]

figure  
imagesc(freq,gate,a');
set(gca,'ColorScale','log')
   colormap summer

  set(gca,'YDir','normal');
  ylabel('Pump Voltage(V)')
  xlabel('Frequency(MHz)')

figure
plot(gate,a1')
  ylabel('Amplitude(arb. units)')
  xlabel('Pump Volatge(V)')
figure
plot(gate,a2')
  ylabel('Amplitude(arb. units)')
  xlabel('Pump Volatge(V)')
