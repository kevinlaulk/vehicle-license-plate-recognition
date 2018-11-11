%test video
clc,clear
close all;
!imaqmem(30000000);               %�����ڴ�ռ�
%=======================��ȡ�豸����==================================
info=imaqhwinfo;%matlab��һ����������
% info =...
%        InstalledAdaptors: {'gentl'  'gige'  'matrox'  'winvideo'}
%         MATLABVersion: '8.1 (R2013a)'
%           ToolboxName: 'Image Acquisition Toolbox'
%        ToolboxVersion: '4.5 (R2013a)'
win_info=imaqhwinfo('winvideo');%winvideo��������������
% win_info =... 
%     AdaptorDllName:'D:\Matlab\toolbox\imaq\imaqadaptors\win64\mwwinvideoimaq.dll' %������dll�ļ�����·��
%     AdaptorDllVersion: '4.5 (R2013a)'%������dll�ļ��汾
%           AdaptorName: 'winvideo'%����������
%             DeviceIDs: {[1]}%�豸ID��
%            DeviceInfo: [1x1 struct]%�豸��Ϣ
win_info.DeviceIDs;%�豸ID��
dev_win_info=win_info.DeviceInfo;%�豸��Ϣ
% dev_win_info =... 
%        DefaultFormat: 'MJPG_1280x720'%��ȡͼƬ��Ĭ�ϸ�ʽ
%        DeviceFileSupported: 0
%                 DeviceName: 'Lenovo EasyCamera'%�豸����
%                   DeviceID: 1%�豸��
%      VideoInputConstructor: 'videoinput('winvideo', 1)'%���󹹽���ʽ��������󲿷ֶ���һ����
%     VideoDeviceConstructor: 'imaq.VideoDevice('winvideo', 1)'
%           SupportedFormats: {1x10 cell}%��ȡ��ͼ��֧�ָ�ʽ
%============================��ȡ����=================================

%============================��ƵԤ���ɼ�����========================
% win_adaptorname=win_info.AdaptorName;%��������������  'winvideo'
% win_deviceID=dev_win_info.DeviceID;%�����������豸��   1
% win_format='MJPG_640x480';%SupportedFormats��ʽ(4)Ϊ  MJPG_640x480
obj = videoinput('winvideo',2,'MJPG_640x480');%������Ƶ�������
obj = videoinput('winvideo',1,'YUY2_640x480');%������Ƶ�������
triggerconfig(obj,'manual');%��Ϊ������ʹ�ܹ�ʵʱ��ʾ
start(obj);
%����ͷԤ������
%web -browser http://blog.csdn.net/msq19895070/article/details/8019663
% vidRes = get(obj, 'VideoResolution'); %��ȡ����
% nBands = get(obj, 'NumberOfBands'); %��ȡ����
% figure()%ָ��Ԥ��������ʾ��figure
% axes()%ָ��Ԥ��������ʾ������ϵ
% hImage = image( zeros(vidRes(2), vidRes(1), nBands) ); 
% preview(obj, hImage);%����ƵԤ������

%����ͷ��ʾ���ó���
% h=figure('NumberTitle','off','Name','��Ƶ',...
%     'MenuBar','none','color','c',...
%     'Position', [0, 0, 1, 1], 'Visible', 'on');         %�½�����
% set(h,'doublebuffer','on','outerposition',get(0,'screensize'));
% h1=axes('Position', [0.02, 0.1, 0.4, 0.8],'Parent',h); %�½���ʾ����
% hold on;
% axis off;

flag_video=1;%����ͷ�����־
first_flag=1;%���δ����־
count2=0;%��ʱ����
sum_percent=0.2;%������ֵ
sum_percent=sum_percent*640;%��ֵ����
manual_flag=0;%�˹�����

while flag_video==1&&count2<1000
frame = getsnapshot(obj);%��ͼ
frame=ycbcr2rgb(frame);
count2
imshow(frame);hold on;plot([1,640],[240,240]);hold off;title('�Զ��������ճ���');drawnow;%��ʾ��ͼ
if first_flag==1
    frame0=im2bw(frame);%������ֵȷ�� %תΪ��ֵͼ
    sum_frame0=sum(frame0(240,:),2);%��240�ж�ֵͼ���
    sum_frame1=sum_frame0;%��ֹ����ʶ���¼�
    first_flag=2;%�ڶ�֡����
elseif first_flag==2
    frame1 = im2bw(frame);%������ֵȷ�� %תΪ��ֵͼ
    sum_frame1=sum(frame1(240,:),2);
    first_flag=0;%�Ժ�֡����
else
    sum_frame0=sum_frame1;
    frame1 = im2bw(frame);%������ֵȷ�� %תΪ��ֵͼ
    sum_frame1=sum(frame1(240,:),2);
end
sum_diff=abs(sum_frame1-sum_frame0);
if sum_diff>sum_percent||manual_flag==1
    count2=1000
end
count2=count2+1;
end
imwrite(frame,'..\���Ʊ�׼ͼ\snapshot.jpg','jpg');%�����ͼ
stop(obj)
delete(obj);




