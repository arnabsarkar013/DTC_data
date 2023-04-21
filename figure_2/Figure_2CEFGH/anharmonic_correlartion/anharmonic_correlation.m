clear all;
   V_mod=[0:20:360 0]
   ki=30;kf=49;


m=1;
for k=ki:kf  %file numbers over which data has to be collected

s1=strcat('grap_',num2str(k),'\Freq1.csv'); %file path 
s2=strcat('grap_',num2str(k),'\Freq2.csv'); %file path 

data1=dlmread(s1);
data2=dlmread(s2);
quad_X1=data1(:,2);%X quadrature
quad_Y1=data1(:,3);%Y quadrature
quad_X2=data2(:,2);%X quadrature
quad_Y2=data2(:,3);%Y quadrature

data1 = [quad_X1 quad_Y1];
data2 = [quad_X2 quad_Y2];
%------------------------------------------
%Mean of quadratures
   %mean_X=mean(quad_X)./(-7.037e-6);
    %mean_X=20.*log10(mean(quad_X)./(-7.037e-6));
    %------------------------------------------------
    %mean of quadratures...
    mean_X1=mean(quad_X1);
    mean_Y1=mean(quad_Y1);
    mean_X2=mean(quad_X2);
    mean_Y2=mean(quad_Y2);
   %mean_Y=mean(quad_Y)./(6.38e-6);
   %mean_Y=20.*log10(mean(quad_Y)./(6.38e-6));
  %Gain--------------------------------------------- 
 G1=sqrt(mean_X1.^2+mean_Y1.^2)
 G2=sqrt(mean_X2.^2+mean_Y2.^2)
 G12=sqrt(mean_X1.^2+mean_Y2.^2)
 G21=sqrt(mean_X2.^2+mean_Y1.^2)
 
 
 quad_X1=quad_X1-mean_X1;%X quadrature-----centring the quadratures to zero.  
quad_Y1=quad_Y1-mean_Y1;%Y quadrature
quad_X2=quad_X2-mean_X2;%X quadrature
quad_Y2=quad_Y2-mean_Y2;%Y quadrature
 
%  N_G=(G-min(G))./(max(G)-min(G));
%error  corresponding to gain-----------------------
err_X1=std(quad_X1);
err_Y1=std(quad_Y1);
err_G1=sqrt(err_X1.^2+err_Y1.^2);
err_X2=std(quad_X2);
err_Y2=std(quad_Y2);
err_G2=sqrt(err_X2.^2+err_Y2.^2);
err_G12=sqrt(err_X1.^2+err_Y2.^2);

err_G21=sqrt(err_X2.^2+err_Y1.^2);

%  N_err_G=(err_G-min(err_G))./(max(err_G)-min(err_G));
%
m_X_t1(m)=sum(mean_X1);
m_Y_t1(m)=sum(mean_Y1);
G1_(m)=sum(G1);
m_X_t2(m)=sum(mean_X2);
m_Y_t2(m)=sum(mean_Y2);
G2_(m)=sum(G2);
G12_(m)=sum(G12);
G21_(m)=sum(G21);
% N_G1(k-ki+1)=sum(N_G);
err_G1_(m)=sum(err_G1);
err_G2_(m)=sum(err_G2);
err_G12_(m)=sum(err_G12);
err_G21_(m)=sum(err_G21);
% N_err_G1(k-ki+1)=sum(N_err_G);
 angle_D1=atan2(mean_Y1,mean_X1)
angle_D1=angle_D1.*(180./pi);
angle_D1_(m)= sum(angle_D1)

 angle_D2=atan2(mean_Y2,mean_X2)
angle_D2=angle_D2.*(180./pi);
angle_D2_(m)= sum(angle_D2)

 angle_D12=atan2(mean_Y2,mean_X1)
angle_D12=angle_D12.*(180./pi);
angle_D12_(m)=sum(angle_D12);

angle_D21=atan2(mean_Y1,mean_X2)
angle_D21=angle_D21.*(180./pi);
angle_D21_(m)=sum(angle_D21);
%Normalisation(b/w 0 to 1)
% N_X=(quad_X-min(quad_X))./(max(quad_X)-min(quad_X));
% N_Y=(quad_Y-min(quad_Y))./(max(quad_Y)-min(quad_Y));
N_X1=quad_X1;
N_Y1=quad_Y1;
N_X2=quad_X2;
N_Y2=quad_Y2;


Nm_X_t1(m)=sum(N_X1);
Nm_Y_t1(m)=sum(N_Y1);
Nm_X_t2(m)=sum(N_X2);
Nm_Y_t2(m)=sum(N_Y2);
% std_NX=std(N_X)% standard deviation
% std_NY=std(N_Y)
% std_NX=std(N_X)% standard deviation
%  std_NY=std(N_Y)
%shifting to zero
 N_X1=N_X1;
N_Y1=N_Y1;
data1 = [N_X1 N_Y1];
 N_X2=N_X2;
N_Y2=N_Y2;
%  N_X1=N_X1-mean(N_X1);
% N_Y1=N_Y1-mean(N_Y1);
% data1 = [N_X1 N_Y1];
%  N_X2=N_X2-mean(N_X2);
% N_Y2=N_Y2-mean(N_Y2);
data2 = [N_X2 N_Y2];
data12 = [N_X1 N_Y2];
data21 = [N_X2 N_Y1];
%-------------------------------------------
% Calculate the eigenvectors and eigenvalues
covariance1 = cov(data1);
covariance2 = cov(data2);
covariance12 = cov(data12);
covariance21 = cov(data21);
[eigenvec_1, eigenval_1 ] = eig(covariance1);
[eigenvec_2, eigenval_2 ] = eig(covariance2);
[eigenvec_12, eigenval_12 ] = eig(covariance12);
[eigenvec_21, eigenval_21 ] = eig(covariance21);
% Get the index of the largest eigenvector
[largest_eigenvec_ind_c_1, r_1] = find(eigenval_1 == max(max(eigenval_1)));
largest_eigenvec_1 = eigenvec_1(:, largest_eigenvec_ind_c_1);
[largest_eigenvec_ind_c_2, r_2] = find(eigenval_2 == max(max(eigenval_2)));
largest_eigenvec_2 = eigenvec_2(:, largest_eigenvec_ind_c_2);
[largest_eigenvec_ind_c_12, r_12] = find(eigenval_12 == max(max(eigenval_12)));
largest_eigenvec_12 = eigenvec_12(:, largest_eigenvec_ind_c_12);
[largest_eigenvec_ind_c_21, r_21] = find(eigenval_21 == max(max(eigenval_21)));
largest_eigenvec_21 = eigenvec_21(:, largest_eigenvec_ind_c_21);
% Get the largest eigenvalue
largest_eigenval_1 = max(max(eigenval_1));
largest_eigenval_2 = max(max(eigenval_2));
largest_eigenval_12 = max(max(eigenval_12));
largest_eigenval_21 = max(max(eigenval_21));
% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c_1 == 1)
    smallest_eigenval_1 = max(eigenval_1(:,2));
    smallest_eigenvec_1 = eigenvec_1(:,2);
else
    smallest_eigenval_1 = max(eigenval_1(:,1));
    smallest_eigenvec_1 = eigenvec_1(1,:);
end
if(largest_eigenvec_ind_c_2 == 1)
    smallest_eigenval_2 = max(eigenval_2(:,2));
    smallest_eigenvec_2 = eigenvec_2(:,2);
else
    smallest_eigenval_2 = max(eigenval_2(:,1));
    smallest_eigenvec_2 = eigenvec_2(1,:);
end
if(largest_eigenvec_ind_c_12 == 1)
    smallest_eigenval_12 = max(eigenval_12(:,2));
    smallest_eigenvec_12 = eigenvec_12(:,2);
else
    smallest_eigenval_12 = max(eigenval_12(:,1));
    smallest_eigenvec_12 = eigenvec_12(1,:);
end
if(largest_eigenvec_ind_c_21 == 1)
    smallest_eigenval_21 = max(eigenval_21(:,2));
    smallest_eigenvec_21 = eigenvec_21(:,2);
else
    smallest_eigenval_21 = max(eigenval_21(:,1));
    smallest_eigenvec_21 = eigenvec_21(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle_1 = atan2(largest_eigenvec_1(2), largest_eigenvec_1(1));
angle_2 = atan2(largest_eigenvec_2(2), largest_eigenvec_2(1));
angle_12 = atan2(largest_eigenvec_12(2), largest_eigenvec_12(1));
angle_21 = atan2(largest_eigenvec_21(2), largest_eigenvec_21(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle_1 < 0)
    angle_1 = angle_1 + 2*pi;
end
% angle_1=angle_1.*(180./pi)
if(angle_2 < 0)
    angle_2 = angle_2 + 2*pi;
end
% angle_2=angle_2.*(180./pi)
if(angle_12 < 0)
    angle_12 = angle_12 + 2*pi;
end
% angle_12=angle_12.*(180./pi)
if(angle_21 < 0)
    angle_21 = angle_21 + 2*pi;
end
% angle_21=angle_21.*(180./pi)

% Get the coordinates of the data mean
avg1 = mean(data1);
avg2 = mean(data2);
avg12 = mean(data12);
avg21 = mean(data21);

% Get the 95% confidence interval error ellipse
% chisquare_val = sqrt(5.991).;
chisquare_val =1.;
theta_grid = linspace(0,2.*pi);
phi_1 = angle_1;
phi_2 = angle_2;
phi_12 = angle_12;
phi_21 = angle_21;
X0_1=avg1(1);
Y0_1=avg1(2);
X0_2=avg2(1);
Y0_2=avg2(2);
X0_12=avg12(1);
Y0_12=avg12(2);
X0_21=avg21(1);
Y0_21=avg21(2);
a1=chisquare_val.*sqrt(largest_eigenval_1);
b1=chisquare_val.*sqrt(smallest_eigenval_1);
a2=chisquare_val.*sqrt(largest_eigenval_2);
b2=chisquare_val.*sqrt(smallest_eigenval_2);
a12=chisquare_val.*sqrt(largest_eigenval_12);
b12=chisquare_val.*sqrt(smallest_eigenval_12);
a21=chisquare_val.*sqrt(largest_eigenval_21);
b21=chisquare_val.*sqrt(smallest_eigenval_21);
% the ellipse in x and y coordinates
ellipse_x_r_1  = a1.*cos( theta_grid );
ellipse_y_r_1  = b1.*sin( theta_grid );
ellipse_x_r_2  = a2.*cos( theta_grid );
ellipse_y_r_2  = b2.*sin( theta_grid );
ellipse_x_r_12  = a12.*cos( theta_grid );
ellipse_y_r_12  = b12.*sin( theta_grid );
ellipse_x_r_21  = a21.*cos( theta_grid );
ellipse_y_r_21  = b21.*sin( theta_grid );
%Define a rotation matrix
R_1 = [cos(phi_1) sin(phi_1); -sin(phi_1) cos(phi_1)];
R_2 = [ cos(phi_2) sin(phi_2); -sin(phi_2) cos(phi_2)];
R_12 = [ cos(phi_12) sin(phi_12); -sin(phi_12) cos(phi_12)];
R_21 = [ cos(phi_21) sin(phi_21); -sin(phi_21) cos(phi_21)];


%let's rotate the ellipse to some angle phi
r_ellipse_1 = [ellipse_x_r_1;ellipse_y_r_1]' *R_1;
r_ellipse_2 = [ellipse_x_r_2;ellipse_y_r_2]' *R_2;
r_ellipse_12 = [ellipse_x_r_12;ellipse_y_r_12]'*R_12;
r_ellipse_21 = [ellipse_x_r_21;ellipse_y_r_21]'*R_21;


figure
subplot(2,2,1)
% Draw the error ellipse
plot(r_ellipse_1(:,1) + X0_1,r_ellipse_1(:,2) + Y0_1,'r','LineWidth',2)
hold on;
% Plot the original data
plot(data1(:,1), data1(:,2), 'r.');
mindata_1 = min(min(data1));
maxdata_1 = max(max(data1));
% Plot the eigenvectors
quiver(X0_1, Y0_1, largest_eigenvec_1(1).*sqrt(largest_eigenval_1), largest_eigenvec_1(2).*sqrt(largest_eigenval_1), '-m', 'LineWidth',2);
quiver(X0_1, Y0_1, smallest_eigenvec_1(1).*sqrt(smallest_eigenval_1), smallest_eigenvec_1(2).*sqrt(smallest_eigenval_1), '-k', 'LineWidth',2);
% hold on;
% Set the axis labels
hXLabel = xlabel('X1');
hYLabel = ylabel('Y1');
hold off;

subplot(2,2,2)
% Draw the error ellipse
plot(r_ellipse_2(:,1) + X0_2,r_ellipse_2(:,2) + Y0_2,'r','LineWidth',2)
hold on;
% Plot the original data
plot(data2(:,1), data2(:,2), 'k.');
mindata_2 = min(min(data2));
maxdata_2 = max(max(data2));
% Plot the eigenvectors
quiver(X0_2, Y0_2, largest_eigenvec_2(1).*sqrt(largest_eigenval_2), largest_eigenvec_2(2).*sqrt(largest_eigenval_2), '-m', 'LineWidth',2);
quiver(X0_2, Y0_2, smallest_eigenvec_2(1).*sqrt(smallest_eigenval_2), smallest_eigenvec_2(2).*sqrt(smallest_eigenval_2), '-k', 'LineWidth',2);
% hold on;
% Set the axis labels
hXLabel = xlabel('X2');
hYLabel = ylabel('Y2');
hold off;

subplot(2,2,3)
% Draw the error ellipse
plot(r_ellipse_12(:,1) + X0_12,r_ellipse_12(:,2) + Y0_12,'r','LineWidth',2)
hold on;
% Plot the original data
plot(data12(:,1), data12(:,2), 'g.');
mindata_12 = min(min(data12));
maxdata_12 = max(max(data12));
% Plot the eigenvectors
quiver(X0_12, Y0_12, largest_eigenvec_12(1).*sqrt(largest_eigenval_12), largest_eigenvec_12(2).*sqrt(largest_eigenval_12), '-m', 'LineWidth',2);
quiver(X0_12, Y0_12, smallest_eigenvec_12(1).*sqrt(smallest_eigenval_12), smallest_eigenvec_12(2).*sqrt(smallest_eigenval_12), '-k', 'LineWidth',2);
% hold on;
% Set the axis labels
hXLabel = xlabel('X_1');
hYLabel = ylabel('Y_2');
hold off;

subplot(2,2,4)
% Draw the error ellipse
plot(r_ellipse_21(:,1) + X0_21,r_ellipse_21(:,2) + Y0_21,'r','LineWidth',2)
hold on;
% Plot the original data
plot(data21(:,1), data21(:,2), 'k.');
mindata_21 = min(min(data21));
maxdata_21 = max(max(data21));
% Plot the eigenvectors
quiver(X0_21, Y0_21, largest_eigenvec_21(1).*sqrt(largest_eigenval_21), largest_eigenvec_21(2).*sqrt(largest_eigenval_21), '-m', 'LineWidth',2);
quiver(X0_21, Y0_21, smallest_eigenvec_21(1).*sqrt(smallest_eigenval_21), smallest_eigenvec_21(2).*sqrt(smallest_eigenval_21), '-k', 'LineWidth',2);
% hold on;
% Set the axis labels
hXLabel = xlabel('X_2');
hYLabel = ylabel('Y_1');
hold off;

%std of mode 1
a1_t(m)=sum(a1);
b1_t(m)=sum(b1);

%std of mode 2
a2_t(m)=sum(a2);
b2_t(m)=sum(b2);

a12_t(m)=sum(a12);
b12_t(m)=sum(b12);
a21_t(m)=sum(a21);
b21_t(m)=sum(b21);
angle1_t(m)=sum(angle_1);
angle2_t(m)=sum(angle_2);
angle12_t(m)=sum(angle_12);
angle21_t(m)=sum(angle_21);
m=m+1;
end
% a_tp=a_t';
% b_tp=b_t';
G1p=G1_';
G2p=G2_';
G12p=G12_';
G21p=G21_';
err_G1p=err_G1_';
err_G2p=err_G2_';
    