function varargout = videoGUI(varargin)
% VIDEOGUI MATLAB code for videoGUI.fig
%      VIDEOGUI, by itself, creates a new VIDEOGUI or raises the existing
%      singleton*.
%
%      H = VIDEOGUI returns the handle to a new VIDEOGUI or the handle to
%      the existing singleton*.
%
%      VIDEOGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIDEOGUI.M with the given input arguments.
%
%      VIDEOGUI('Property','Value',...) creates a new VIDEOGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before videoGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to videoGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help videoGUI

% Last Modified by GUIDE v2.5 17-Mar-2016 21:58:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @videoGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @videoGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before videoGUI is made visible.
function videoGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to videoGUI (see VARARGIN)

% Choose default command line output for videoGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes videoGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = videoGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%选择视频
%设置默认值
global str
global flag_video
global count_value
global sum_percent1%触发阈值
global manual_flag%人工触发
global erzhi_yuzhi%二值化阈值方法
global tezheng
global shibiefangfa
global zongshibie %总识别
flag_video=1;%摄像头允许标志
sum_percent1=0.05;%触发阈值
manual_flag=0;%关闭人工触发
count_value=100000;%默认不计时触发，使为无穷大
tezheng=2;%结构
shibiefangfa=1;%模板匹配
erzhi_yuzhi=3;%综合阈值
zongshibie =1;%识别方法（1080*1920）
[filename,pathname]=...
    uigetfile({'*.avi';'*.mp4';'*.flv' },'选择视频');
str=[pathname filename];

% info=aviinfo(str);
% fum=info.NumFrames;
% count=0;%计数
% for i=1:fum
%     mov=aviread(str,i);
%     I=mov.cdata;
%     J=rgb2gray(I);
%     count=count+1
% end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%停止监控
global flag_video
flag_video=0;
clc,clear;

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%退出
clc,clear;
close(gcf)

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
%定时识别
global count_value
str=get(handles.popupmenu1,'value');
switch str
    case 1
        count_value=1000;
    case 2
        count_value=100;
    case 3
        count_value=500;
    case 4
        count_value=2000;
    case 5
        count_value=5000;
    case 6
        count_value=10000;
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
global sum_percent1
str=get(handles.popupmenu2,'value');
switch str
    case 1
        sum_percent1=0.01;
    case 2
        sum_percent1=0.02;
    case 3
        sum_percent1=0.03;
    case 4
        sum_percent1=0.05;
    case 5
        sum_percent1=0.08;
    case 6
        sum_percent1=0.10;
    case 7
        sum_percent1=0.15;
    case 8
        sum_percent1=0.20;
    case 9
        sum_percent1=0.25;
    case 10
        sum_percent1=0.50;
    case 11
        sum_percent1=0.80;
end


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%人工触发
global manual_flag
manual_flag=1;

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%播放视频
global str
global flag_video%是否结束视频
global count_value%计数触发
global sum_percent1%触发阈值
global manual_flag%人工触发
global im
global erzhi_yuzhi%二值化阈值方法
global tezheng    %bp网络识别方法
global shibiefangfa %识别方法
global zongshibie%总识别方法
global RegCode
height=1080;%视频高度
height1=height/2;
width=1920;%视频宽度
first_flag=1;%初次储存标志
sum_percent=sum_percent1*width;%阈值计算
count2=0;%计时触发识别
mov = VideoReader( str );%读取视频
for i=1:mov.NumberOfFrames
    frame = read( mov, i );
    count2 %显示触发时间
    axes(handles.axes5);imshow(frame);hold on;plot([1,width],[height1,height1]);hold off;drawnow;%显示截图
if first_flag==1
    frame0=im2bw(frame);%触发阈值确定 %转为二值图
    sum_frame0=sum(frame0(height1,:),2);%触发行二值图求和
    sum_frame1=sum_frame0;%防止触发识别事件
    first_flag=2;%第二帧处理
elseif first_flag==2
    frame1 = im2bw(frame);%触发阈值确定 %转为二值图
    sum_frame1=sum(frame1(height1,:),2);
    first_flag=0;%以后帧处理
else
    sum_frame0=sum_frame1;
    frame1 = im2bw(frame);%触发阈值确定 %转为二值图
    sum_frame1=sum(frame1(height1,:),2);
    axes(handles.axes6);imshow(frame1);hold on;plot([1,width],[height1,height1]);hold off;drawnow;%显示二值截图
end

sum_diff=abs(sum_frame1-sum_frame0)
if sum_diff>sum_percent||manual_flag==1||count2==count_value
    manual_flag=0;%关闭人工触发
    im=frame;
%     [M,N,C]=size(im);%检查读入帧格式
    axes(handles.axes7);imshow(im);%执行函数
    file_name1=strcat('..\车牌标准图\snapshot','.jpg');%截图路径
    imwrite(frame,file_name1,'jpg');%保存截图
    if zongshibie==1
        shibiechengxu2();
    elseif zongshibie==2
        blue_recog();%执行车牌辨别
    else
        shibiechengxu3();
    end
    %======显示结果=============
for l=1:7
file_name1=strcat('..\切割子图\切割子图',int2str(l),'.jpg');
SegGray1=imread(file_name1);
switch l
    case 1
        axes(handles.axes8);
    case 2
        axes(handles.axes9);
    case 3
        axes(handles.axes10);
    case 4
        axes(handles.axes11);
    case 5
        axes(handles.axes12);
    case 6 
        axes(handles.axes13);
    case 7
        axes(handles.axes14);
end
imshow(SegGray1);
end
      %======显示字符===========
set(handles.edit1,'string',RegCode);

end
count2=count2+1;
if flag_video==0
    return;
end
end
clc,clear

% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global erzhi_yuzhi
str=get(hObject,'string');
switch str
    case '平均最佳阈值'
        erzhi_yuzhi=1;
    case 'OSTU'
        erzhi_yuzhi=2;
    case '综合最佳阈值'
        erzhi_yuzhi=3;
end


% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
%识别模式选择

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Menu 作者菜单
msgbox('programer：kevin  ','Author','')

% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Menu 帮助菜单
help_text1();

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%打开文件
global str
[filename,pathname]=...
    uigetfile({'*.avi';'*.mp4';'*.flv'},'选择视频');
str=[pathname filename];

% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%播放
pushbutton5_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%人工触发
global manual_flag
manual_flag=1;

% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%退出
clc,clear;
close(gcf)


% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global shibiefangfa
str=get(hObject,'string');
switch str
    case '模板识别'
        shibiefangfa=1;
    case 'BP神经网络'
        shibiefangfa=2;
    case '欧氏距离'
        shibiefangfa=3;
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel4.
function uipanel4_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel4 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global tezheng 
str=get(hObject,'string');
switch str
    case '统计特征'
        tezheng=1;
    case '结构特征'
        tezheng=2;
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
global zongshibie
str=get(handles.popupmenu3,'value');
switch str
    case 1
        zongshibie=1
    case 2
        zongshibie=2
    case 3
        zongshibie=3
end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
