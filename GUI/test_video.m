%test video
clc,clear
close all;
!imaqmem(30000000);               %申请内存空间
%=======================获取设备参数==================================
info=imaqhwinfo;%matlab第一适配器参数
% info =...
%        InstalledAdaptors: {'gentl'  'gige'  'matrox'  'winvideo'}
%         MATLABVersion: '8.1 (R2013a)'
%           ToolboxName: 'Image Acquisition Toolbox'
%        ToolboxVersion: '4.5 (R2013a)'
win_info=imaqhwinfo('winvideo');%winvideo第四适配器参数
% win_info =... 
%     AdaptorDllName:'D:\Matlab\toolbox\imaq\imaqadaptors\win64\mwwinvideoimaq.dll' %适配器dll文件绝对路径
%     AdaptorDllVersion: '4.5 (R2013a)'%适配器dll文件版本
%           AdaptorName: 'winvideo'%适配器名称
%             DeviceIDs: {[1]}%设备ID号
%            DeviceInfo: [1x1 struct]%设备信息
win_info.DeviceIDs;%设备ID号
dev_win_info=win_info.DeviceInfo;%设备信息
% dev_win_info =... 
%        DefaultFormat: 'MJPG_1280x720'%获取图片的默认格式
%        DeviceFileSupported: 0
%                 DeviceName: 'Lenovo EasyCamera'%设备名称
%                   DeviceID: 1%设备号
%      VideoInputConstructor: 'videoinput('winvideo', 1)'%对象构建方式，这个绝大部分都是一样的
%     VideoDeviceConstructor: 'imaq.VideoDevice('winvideo', 1)'
%           SupportedFormats: {1x10 cell}%获取的图像支持格式
%============================获取结束=================================

%============================视频预览采集保存========================
% win_adaptorname=win_info.AdaptorName;%第四适配器名称  'winvideo'
% win_deviceID=dev_win_info.DeviceID;%第四适配器设备号   1
% win_format='MJPG_640x480';%SupportedFormats格式(4)为  MJPG_640x480
obj = videoinput('winvideo',2,'MJPG_640x480');%创建视频输入对象
obj = videoinput('winvideo',1,'YUY2_640x480');%创建视频输入对象
triggerconfig(obj,'manual');%人为触发，使能够实时显示
start(obj);
%摄像头预览程序
%web -browser http://blog.csdn.net/msq19895070/article/details/8019663
% vidRes = get(obj, 'VideoResolution'); %获取参数
% nBands = get(obj, 'NumberOfBands'); %获取参数
% figure()%指定预览窗体显示的figure
% axes()%指定预览窗口显示的坐标系
% hImage = image( zeros(vidRes(2), vidRes(1), nBands) ); 
% preview(obj, hImage);%打开视频预览窗口

%摄像头显示设置程序
% h=figure('NumberTitle','off','Name','视频',...
%     'MenuBar','none','color','c',...
%     'Position', [0, 0, 1, 1], 'Visible', 'on');         %新建窗口
% set(h,'doublebuffer','on','outerposition',get(0,'screensize'));
% h1=axes('Position', [0.02, 0.1, 0.4, 0.8],'Parent',h); %新建显示窗口
% hold on;
% axis off;

flag_video=1;%摄像头允许标志
first_flag=1;%初次储存标志
count2=0;%计时触发
sum_percent=0.2;%触发阈值
sum_percent=sum_percent*640;%阈值计算
manual_flag=0;%人工触发

while flag_video==1&&count2<1000
frame = getsnapshot(obj);%截图
frame=ycbcr2rgb(frame);
count2
imshow(frame);hold on;plot([1,640],[240,240]);hold off;title('自动触发拍照程序');drawnow;%显示截图
if first_flag==1
    frame0=im2bw(frame);%触发阈值确定 %转为二值图
    sum_frame0=sum(frame0(240,:),2);%第240行二值图求和
    sum_frame1=sum_frame0;%防止触发识别事件
    first_flag=2;%第二帧处理
elseif first_flag==2
    frame1 = im2bw(frame);%触发阈值确定 %转为二值图
    sum_frame1=sum(frame1(240,:),2);
    first_flag=0;%以后帧处理
else
    sum_frame0=sum_frame1;
    frame1 = im2bw(frame);%触发阈值确定 %转为二值图
    sum_frame1=sum(frame1(240,:),2);
end
sum_diff=abs(sum_frame1-sum_frame0);
if sum_diff>sum_percent||manual_flag==1
    count2=1000
end
count2=count2+1;
end
imwrite(frame,'..\车牌标准图\snapshot.jpg','jpg');%保存截图
stop(obj)
delete(obj);




