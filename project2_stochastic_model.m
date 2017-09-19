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

novatel_wz_bias = -2.3/(180*3600/pi); % converted deg/hr to rad/sec;
novatel_wz_sf = -20000/1000000; %converted ppm to actual SF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xbow_fx_ss = (xbow_fx - xbow_fx_bias)/(1+xbow_fx_sf);  %%%%%%%%%%  RESIDUAL??
xbow_fy_ss = (xbow_fy - xbow_fy_bias)/(1+xbow_fy_sf);
xbow_fz_ss = (xbow_fz - xbow_fz_bias)/(1+xbow_fz_sf);


xbow_wx_ss = (xbow_wx - xbow_wx_bias)/(1+xbow_wx_sf);  %%%%%%%%%%  RESIDUAL??
xbow_wy_ss = (xbow_wy - xbow_wy_bias)/(1+xbow_wy_sf);
xbow_wz_ss = (xbow_wz - xbow_wz_bias)/(1+xbow_wz_sf);


%ac_xbow_fx_ss = xcorr2(xbow_fx_ss(1:3000000));
%ac_xbow_fy_ss = xcorr2(xbow_fy_ss(1:3000000));
%ac_xbow_fz_ss = xcorr2(xbow_fz_ss(1:3000000));

%ac_xbow_wx_ss = xcorr2(xbow_wx_ss(1:3000000));
%ac_xbow_wy_ss = xcorr2(xbow_wy_ss(1:3000000));
%ac_xbow_wz_ss = xcorr2(xbow_wz_ss(1:3000000));

xbow_wz_ss = xbow_fz_ss;
span = 50;
smooth_xbow_wz_ss_50 = smooth(xbow_wz_ss,span);
std_smooth(1) = std(smooth_xbow_wz_ss_50);

span = 100;
smooth_xbow_wz_ss_100 = smooth(xbow_wz_ss,span);
std_smooth(2) = std(smooth_xbow_wz_ss_100);

span = 200;
smooth_xbow_wz_ss_200 = smooth(xbow_wz_ss,span);
std_smooth(3) = std(smooth_xbow_wz_ss_200);

span = 400;
smooth_xbow_wz_ss_400 = smooth(xbow_wz_ss,span);
std_smooth(4) = std(smooth_xbow_wz_ss_400);

span = 600;
smooth_xbow_wz_ss_600 = smooth(xbow_wz_ss,span);
std_smooth(5) = std(smooth_xbow_wz_ss_600);

span = 1000;
smooth_xbow_wz_ss_1000 = smooth(xbow_wz_ss,span);
std_smooth(6) = std(smooth_xbow_wz_ss_1000);

span = 1600;
smooth_xbow_wz_ss_1600 = smooth(xbow_wz_ss,span);
std_smooth(7) = std(smooth_xbow_wz_ss_1600);

span = 2000;
smooth_xbow_wz_ss_2000 = smooth(xbow_wz_ss,span);
std_smooth(8) = std(smooth_xbow_wz_ss_2000);




x_axis = (1/200)*(1:length(xbow_wz_ss));

subplot(2,1,1)
plot(x_axis,smooth_xbow_wz_ss_50,'b'); hold on; plot(x_axis,smooth_xbow_wz_ss_200, 'r'); plot(x_axis,smooth_xbow_wz_ss_400, 'g');hold on; plot(x_axis,smooth_xbow_wz_ss_1000,'y');
     lg=legend('Span=0.25 s','Span=1 s','Span=2 s','Span=5 s');
    gt1=findobj(lg,'type','text');
    set(gt1,'fontname','--','fontweight','bold');
    
    xlabel('time (seconds)','fontweight','bold','fontsize',12);
    ylabel('Xbow Acc-Z Axis Data (m/s^2)','fontweight','bold','fontsize',12);
    title('Xbow Acc denoising with MA','fontweight','bold','fontsize',16)
    
    x_a = [0.25 0.5 1 2 3 5 8 10];
    subplot(2,1,2)
    plot(x_a, std_smooth,'r*-','LineWidth',2);
    xlabel('MA span (seconds)','fontweight','bold','fontsize',12);
    ylabel('Standard deviation - Xbow Acc Z axis (m/s^2)','fontweight','bold','fontsize',12);
    title('Standad deviation variation v/s MA span','fontweight','bold','fontsize',16)






