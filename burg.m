% Ranjeeth KS, University of Calgary
clc;
clear all;
close all;
format long g;
% test location
phi = 44.230697181686281*pi/180;
lambda = -76.466848939548527*pi/180;
height = 95.177099999999996;
gravity	 = 9.805209209982110;

%%%%% GRAVITY MODEL %%%%%
a1=9.7803267714; 
a4=-0.0000030876910891;
a2=0.0052790414; 
a5=0.0000000043977311;
a3=0.0000232718; 
a6=0.0000000000007211;
g = a1*(1+a2*sin(phi)*sin(phi)+a3*sin(phi)*sin(phi)*sin(phi)*sin(phi))+(a4+a5*sin(phi)*sin(phi))*height+a6*height*height;

we = 7.2921151467e-5; % rad/sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long g;
%novatel_data = load('D:\Study\ENGO 623\Project 1 and 2\Project2\static_Nova_16Aug07.mat');
novatel_data = load('static_Nova_16Aug07.mat');

xbow = load('xbow_fx.mat');
xbow_fx = xbow.xbow_fx;

xbow = load('xbow_fy.mat');
xbow_fy = xbow.xbow_fy;

xbow = load('xbow_fz.mat');
xbow_fz = xbow.xbow_fz;

xbow = load('xbow_wx.mat');
xbow_wx = xbow.xbow_wx;

xbow = load('xbow_wy.mat');
xbow_wy = xbow.xbow_wy;

xbow = load('xbow_wz.mat');
xbow_wz = xbow.xbow_wz;

novatel_fx = novatel_data.f.x;
novatel_fy = novatel_data.f.y;
novatel_fz = novatel_data.f.z;

novatel_wx = novatel_data.w.x;
novatel_wy = novatel_data.w.y;
novatel_wz = novatel_data.w.z;

xbow_fx_bias =  4.6*1e-3; %converted mg to g 
xbow_fx_sf = -0.02/100; %converted % to actual SF

xbow_fy_bias =  7.25*1e-3; %converted mg to g 
xbow_fy_sf = -0.066/100; %converted % to actual SF

xbow_fz_bias =  -11.65*1e-3; %converted mg to g 
xbow_fz_sf = -0.025/100; %converted % to actual SF

xbow_wx_bias = 0.079; % No Conversion is needed, already in deg/ sec;
xbow_wx_sf = -0.05/100; %converted % to actual SF

xbow_wy_bias = 0.088; % No Conversion is needed, already in deg/ sec;
xbow_wy_sf = 0.17/100; %converted % to actual SF

xbow_wz_bias = 0.198; % No Conversion is needed, already in deg/ sec;
xbow_wz_sf = 0.17/100; %converted % to actual SF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
novatel_fx_bias =  -0.48*1e-3*g; %converted mg to m/s^2 
novatel_fx_sf = 65/1000000; %converted ppm to actual SF


novatel_fy_bias =  0.165*1e-3*g; %converted mg to m/s^2 
novatel_fy_sf = 45/1000000; %converted ppm to actual SF


novatel_wz_bias = -2.3/(180*3600/pi); % converted deg/hr to rad/sec;
novatel_wz_sf = -20000/1000000; %converted ppm to actual SF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xbow_fx_ss = (xbow_fx - xbow_fx_bias)/(1+xbow_fx_sf);  %%%%%%%%%%  RESIDUAL??
xbow_fy_ss = (xbow_fy - xbow_fy_bias)/(1+xbow_fy_sf);
xbow_fz_ss = (xbow_fz - xbow_fz_bias)/(1+xbow_fz_sf);


xbow_wx_ss = (xbow_wx - xbow_wx_bias)/(1+xbow_wx_sf);  %%%%%%%%%%  RESIDUAL??
xbow_wy_ss = (xbow_wy - xbow_wy_bias)/(1+xbow_wy_sf);
xbow_wz_ss = (xbow_wz - xbow_wz_bias)/(1+xbow_wz_sf);

novatel_fx_ss = (novatel_fx - novatel_fx_bias)/(1+novatel_fx_sf);  %%%%%%%%%%  RESIDUAL??
novatel_wz_ss = (novatel_wz - novatel_wz_bias)/(1+novatel_wz_sf);
novatel_fy_ss = (novatel_fy - novatel_fy_bias)/(1+novatel_fy_sf);

%ac_xbow_fx_ss = xcorr2(xbow_fx_ss(1:25));
%ac_xbow_fy_ss = xcorr2(xbow_fy_ss(1:25));
%ac_xbow_fz_ss = xcorr2(xbow_fz_ss(1:25));

%ac_xbow_wx_ss = xcorr2(xbow_wx_ss(1:25));
%ac_xbow_wy_ss = xcorr2(xbow_wy_ss(1:25));
corr=smooth(xbow_wz_ss(1:20000),20000);
corr_db3=corr;
%[corr_db3, s, b]=wden(xbow_wz_ss(1:20000), 'rigrsure', 's', 'mln', 0, 'db3');

u_n =corr_db3;
order =  1;
tic;
[a,MSE]=aryule(u_n,order);
time (1,1)= toc;
RMSE(1,1) = MSE;


order =  2;
tic;
[a,MSE]=aryule(u_n,order);
time (1,order)= toc;
RMSE(1,order) = MSE;

order =  3;
tic;
[a,MSE]=aryule(u_n,order);
time (1,order)= toc;
RMSE(1,order) = MSE;


order =  4;
tic;
[a,MSE]=aryule(u_n,order);
time (1,order)= toc;
RMSE(1,order) = MSE;

order =  5;
tic;
[a,MSE]=aryule(u_n,order);
time (1,order)= toc;
RMSE(1,order) = MSE;

order =  6;
tic;
[a,MSE]=aryule(u_n,order);
time (1,order)= toc;
RMSE(1,order) = MSE;

order =  7;
tic;
[a,MSE]=aryule(u_n,order);
time (1,order)= toc;
RMSE(1,order) = MSE;

order =  8;
tic;
[a,MSE]=aryule(u_n,order);
time (1,order)= toc;
RMSE(1,order) = MSE;

order =  9;
tic;
[a,MSE]=aryule(u_n,order);
time (1,order)= toc;
RMSE(1,order) = MSE;

order =  10;
tic;
[a,MSE]=aryule(u_n,order);
time (1,order)= toc;
RMSE(1,order) = MSE;

order =  1;
tic;
[a,MSE]=arcov(u_n,order);
time (2,1)= toc;
RMSE(2,1) = MSE;


order =  2;
tic;
[a,MSE]=arcov(u_n,order);
time (2,order)= toc;
RMSE(2,order) = MSE;

order =  3;
tic;
[a,MSE]=arcov(u_n,order);
time (2,order)= toc;
RMSE(2,order) = MSE;


order =  4;
tic;
[a,MSE]=arcov(u_n,order);
time (2,order)= toc;
RMSE(2,order) = MSE;

order =  5;
tic;
[a,MSE]=arcov(u_n,order);
time (2,order)= toc;
RMSE(2,order) = MSE;

order =  6;
tic;
[a,MSE]=arcov(u_n,order);
time (2,order)= toc;
RMSE(2,order) = MSE;

order =  7;
tic;
[a,MSE]=arcov(u_n,order);
time (2,order)= toc;
RMSE(2,order) = MSE;

order =  8;
tic;
[a,MSE]=arcov(u_n,order);
time (2,order)= toc;
RMSE(2,order) = MSE;

order =  9;
tic;
[a,MSE]=arcov(u_n,order);
time (2,order)= toc;
RMSE(2,order) = MSE;

order =  10;
tic;
[a,MSE]=arcov(u_n,order);
time (2,order)= toc;
RMSE(2,order) = MSE;

order =  1;
tic;
[a,MSE]=arburg(u_n,order);
time (3,order)= toc;
RMSE(3,order) = MSE;


order =  2;
tic;
[a,MSE]=arburg(u_n,order);
time(3,order)= toc;
RMSE(3,order) = MSE;

order =  3;
tic;
[a,MSE]=arburg(u_n,order);
time(3,order)= toc;
RMSE(3,order) = MSE;


order =  4;
tic;
[a,MSE]=arburg(u_n,order);
time(3,order)= toc;
RMSE(3,order) = MSE;

order =  5;
tic;
[a,MSE]=arburg(u_n,order);
time(3,order)= toc;
RMSE(3,order) = MSE;

order =  6;
tic;
[a,MSE]=arburg(u_n,order);
time(3,order)= toc;
RMSE(3,order) = MSE;

order =  7;
tic;
[a,MSE]=arburg(u_n,order);
time(3,order)= toc;
RMSE(3,order) = MSE;

order =  8;
tic;
[a,MSE]=arburg(u_n,order);
time(3,order)= toc;
RMSE(3,order) = MSE;

order =  9;
tic;
[a,MSE]=arburg(u_n,order);
time(3,order)= toc;
RMSE(3,order) = MSE;

order =  10; 
tic;
[a,MSE]=arburg(u_n,order);
time(3,order)= toc;
RMSE(3,order) = MSE;
plot(RMSE(1,1:10),'b*-','LineWidth',2); hold on; plot(RMSE(2,1:10),'r*-','LineWidth',2);plot(RMSE(3,1:10),'g*-','LineWidth',2);
  
% plot(xx_a,time(1,7:15),'b*-','LineWidth',2); hold on; plot(xx_a,time(2,7:15),'r*-','LineWidth',2);plot(xx_a,time(3,7:15),'g*-','LineWidth',2);
lg=legend('Yule walker Method','Covariance Method','Burg Method');
    gt1=findobj(lg,'type','text');
     set(gt1,'fontname','--','fontweight','bold');
     
     xlabel('order','fontweight','bold','fontsize',12);
     ylabel('MSE deg/sec','fontweight','bold','fontsize',12);
    title('Time profile for AR estimation methods','fontweight','bold','fontsize',16)
    h=0;
% a_c = a(2:3);
% for i=20:20000;
% estimate_u_n(1,i)=-a_c(1)*u_n(i-1) + -a_c(2)*u_n(i-2) ;
% end
% 
% SE=0;
% for i=20:20000;
%     SE=SE+(u_n(i)-estimate_u_n(1,i))*(u_n(i)-estimate_u_n(1,i));
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%
% RMSE(1) = sqrt(SE/(20000-20));
% order =  3;
% 
% 
% a=arcov(u_n,order);
% 
% a_c = a(2:4);
% 
% for i=20:20000;
% estimate_u_n(2,i)=-a_c(1)*u_n(i-1) + -a_c(2)*u_n(i-2) + -a_c(3)*u_n(i-3);
% end
% y=1; 
% SE=0;
% for i=20:20000;
%     SE=SE+(u_n(i)-estimate_u_n(2,i))*(u_n(i)-estimate_u_n(2,i));
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%
% RMSE(2) = sqrt(SE/(20000-20));%%%%%%%%%%%%%%%%%%%%%%
% 
% order =  4;
% 
% 
% a=arcov(u_n,order);
% 
% a_c = a(2:5);
% for i=20:20000;
% estimate_u_n(3,i)=-a_c(1)*u_n(i-1) + -a_c(2)*u_n(i-2) + -a_c(3)*u_n(i-3) + -a_c(4)*u_n(i-4);
% end
% 
% SE=0;
% for i=20:20000;
%     SE=SE+(u_n(i)-estimate_u_n(3,i))*(u_n(i)-estimate_u_n(3,i));
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%
% RMSE(3) = sqrt(SE/(20000-20));%%%%%%%%%%%%%%%%%%%%%
% 
% order =  5;
% 
% 
% a=arcov(u_n,order);
% 
% a_c = a(2:6);
% 
% for i=20:20000;
% estimate_u_n(4,i)=-a_c(1)*u_n(i-1) + -a_c(2)*u_n(i-2) + -a_c(3)*u_n(i-3) + -a_c(4)*u_n(i-4) + -a_c(5)*u_n(i-5);
% end
% 
% SE=0;
% for i=20:20000;
%     SE=SE+(u_n(i)-estimate_u_n(4,i))*(u_n(i)-estimate_u_n(4,i));
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%
% RMSE(4) = sqrt(SE/(20000-20));
% %%%%%%%%%%%%%%%%%%%%%
% 
% order =  6;
% 
% 
% a=arcov(u_n,order);
% 
% a_c = a(2:7);
% 
% for i=20:20000;
% estimate_u_n(5,i)=-a_c(1)*u_n(i-1) + -a_c(2)*u_n(i-2) + -a_c(3)*u_n(i-3) + -a_c(4)*u_n(i-4) + -a_c(5)*u_n(i-5) + -a_c(6)*u_n(i-6);
% end
% 
% SE=0;
% for i=20:20000;
%     SE=SE+(u_n(i)-estimate_u_n(5,i))*(u_n(i)-estimate_u_n(5,i));
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%
% RMSE(5) = sqrt(SE/(20000-20));
% %%%%%%%%%%%%%%%%%%%%%
% 
% order =  7;
% 
% 
% a=arcov(u_n,order);
% 
% a_c = a(2:8);
% 
% for i=20:20000;
% estimate_u_n(6,i)=-a_c(1)*u_n(i-1) + -a_c(2)*u_n(i-2) + -a_c(3)*u_n(i-3) + -a_c(4)*u_n(i-4) + -a_c(5)*u_n(i-5) + -a_c(6)*u_n(i-6) + -a_c(7)*u_n(i-7);
% end
% 
% SE=0;
% for i=20:20000;
%     SE=SE+(u_n(i)-estimate_u_n(6,i))*(u_n(i)-estimate_u_n(6,i));
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%
% RMSE(6) = sqrt(SE/(20000-20));
%%%%%%%%%%%%%%%%%%%%

% x_axis=(1:25)/100;
% plot(x_axis,u_n(1:25),'LineWidth',2),hold on; plot(x_axis,estimate_u_n(1,:),'r','LineWidth',2);
% lg=legend('Denoised data','Estimated Data');
%     gt1=findobj(lg,'type','text');
%     set(gt1,'fontname','--','fontweight','bold');
%     
%     xlabel('time (seconds)','fontweight','bold','fontsize',12);
%     ylabel('Novatel fx m/s^2','fontweight','bold','fontsize',12);
%     title('-','fontweight','bold','fontsize',16)