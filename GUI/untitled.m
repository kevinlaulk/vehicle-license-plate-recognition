function varargout = untitled(varargin)

% erzhi_yuzhi=3;
% shibiefangfa=1;
% UNTITLED MATLAB code for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 17-Mar-2016 22:11:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
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


% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
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
%��ȡͼƬ
global im
[filename,pathname]=...
    uigetfile({'*.jpg';'*.bmp';'*.gif'},'ѡ��ͼƬ');
str=[pathname filename];
im=imread(str);
axes(handles.axes1);
imshow(im);
%����Ĭ��ֵ
global tezheng
global erzhi_yuzhi
global shibiefangfa
global shuipingxuanzhuan 
global time_weizhi
tezheng=2;%�ṹ
shibiefangfa=1;%ģ��ƥ��
erzhi_yuzhi=3;%�ۺ���ֵ
shuipingxuanzhuan=1;%hough�任
time_weizhi=1;%λ�ó�ʼ��


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%�ر��˳�
clc,clear
close(gcf)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% h=waitbar(0,'���Եȣ����ڿ�ʼʶ���ļ�...');
% waitbar(counter_h/7,h,'��ȴ���ʶ����...');
% waitbar(1,h,'�����');
% delete(h);
global tezheng
global erzhi_yuzhi
global shibiefangfa
global RegCode
global shuipingxuanzhuan 
global time_weizhi
blue_recog();
time_weizhi=time_weizhi+1;
file_name1=strcat('..\���Ʊ�׼ͼ\�Ҷȳ���','.jpg');%�Ҷȳ���ͼ
subcol=imread(file_name1);
axes(handles.axes2);
imshow(subcol);
file_name1=strcat('..\���Ʊ�׼ͼ\���Ʊ��ͼ','.jpg');%���Ʊ��ͼ
subcol=imread(file_name1);
axes(handles.axes12);
imshow(subcol);
file_name1=strcat('..\���Ʊ�׼ͼ\��ֵ����','.jpg');%�洢��ֵͼ
SegGray1=imread(file_name1);
axes(handles.axes3);
imshow(SegGray1);
for l=1:7
file_name1=strcat('..\�и���ͼ\�и���ͼ',int2str(l),'.jpg');
SegGray1=imread(file_name1);
switch l
    case 1
        axes(handles.axes4);
    case 2
        axes(handles.axes5);
    case 3
        axes(handles.axes6);
    case 4
        axes(handles.axes7);
    case 5
        axes(handles.axes8);
    case 6 
        axes(handles.axes9);
    case 7
        axes(handles.axes10);
end
imshow(SegGray1);
end

set(handles.edit2,'string',RegCode);
clc,clear


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
%Menu ���߲˵�
msgbox('programer��kevin  ','Author','')

% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Menu �����˵�
help_text1();


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Menu �򿪲˵�
global im
[filename,pathname]=...
    uigetfile({'*.jpg';'*.bmp';'*.gif'},'ѡ��ͼƬ');
str=[pathname filename];
im=imread(str);
axes(handles.axes1);
imshow(im);
global tezheng
global erzhi_yuzhi
global shibiefangfa
global shuipingxuanzhuan 
global time_weizhi
tezheng=2;
shibiefangfa=1;
erzhi_yuzhi=3;
shuipingxuanzhuan=1;
time_weizhi=1;

% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Menu ʶ��˵�
global tezheng
global erzhi_yuzhi
global shibiefangfa
global shuipingxuanzhuan 
blue_recog();

file_name1=strcat('..\���Ʊ�׼ͼ\�Ҷȳ���','.jpg');%�洢��׼��ͼ
subcol=imread(file_name1);
axes(handles.axes2);
imshow(subcol);
file_name1=strcat('..\���Ʊ�׼ͼ\��ֵ����','.jpg');%�洢��׼��ͼ
SegGray1=imread(file_name1);
axes(handles.axes3);
imshow(SegGray1);

clc,clear

% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Menu �˳��˵�
clc,clear
close(gcf)


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
    case 'ƽ�������ֵ'
        erzhi_yuzhi=1;
    case 'OSTU'
        erzhi_yuzhi=2;
    case '�ۺ������ֵ'
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
global shibiefangfa
str=get(hObject,'string');
switch str
    case 'ģ��ʶ��'
        shibiefangfa=1;
    case 'BP������'
        shibiefangfa=2;
    case 'ŷ�Ͼ���'
        shibiefangfa=3;
end


% --------------------------------------------------------------------
function uipanel1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tezheng
bp_train();%�������¼

% --- Executes when selected object is changed in uipanel6.
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global tezheng 
str=get(hObject,'string');
switch str
    case 'ͳ������'
        tezheng=1;
    case '�ṹ����'
        tezheng=2;
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel7.
function uipanel7_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel7 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global shuipingxuanzhuan 
str=get(hObject,'string');
switch str
    case 'hough���'
        shuipingxuanzhuan=1;
    case '͸�ӱ任'
        shuipingxuanzhuan=2;
end
