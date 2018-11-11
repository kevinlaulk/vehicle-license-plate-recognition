function varargout = help_text1(varargin)
% HELP_TEXT1 MATLAB code for help_text1.fig
%      HELP_TEXT1, by itself, creates a new HELP_TEXT1 or raises the existing
%      singleton*.
%
%      H = HELP_TEXT1 returns the handle to a new HELP_TEXT1 or the handle to
%      the existing singleton*.
%
%      HELP_TEXT1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HELP_TEXT1.M with the given input arguments.
%
%      HELP_TEXT1('Property','Value',...) creates a new HELP_TEXT1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before help_text1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to help_text1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help help_text1

% Last Modified by GUIDE v2.5 24-Jan-2016 17:46:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @help_text1_OpeningFcn, ...
                   'gui_OutputFcn',  @help_text1_OutputFcn, ...
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


% --- Executes just before help_text1 is made visible.
function help_text1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to help_text1 (see VARARGIN)

% Choose default command line output for help_text1
handles.output = hObject;

%自动打开帮助文件
fin=fopen('project_recode.txt','r');
f=fscanf(fin,'%s');
set(handles.edit1,'string',f);
fclose(fin);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes help_text1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = help_text1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%关闭该GUI
clc,clear
close(gcf);
