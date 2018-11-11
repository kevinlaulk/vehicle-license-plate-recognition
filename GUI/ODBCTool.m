function varargout = ODBCTool(varargin)
% ODBCTOOL MATLAB code for ODBCTool.fig
%      ODBCTOOL, by itself, creates a new ODBCTOOL or raises the existing
%      singleton*.
%
%      H = ODBCTOOL returns the handle to a new ODBCTOOL or the handle to
%      the existing singleton*.
%
%      ODBCTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ODBCTOOL.M with the given input arguments.
%
%      ODBCTOOL('Property','Value',...) creates a new ODBCTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ODBCTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ODBCTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ODBCTool

% Last Modified by GUIDE v2.5 22-Mar-2016 08:14:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ODBCTool_OpeningFcn, ...
                   'gui_OutputFcn',  @ODBCTool_OutputFcn, ...
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


% --- Executes just before ODBCTool is made visible.
function ODBCTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ODBCTool (see VARARGIN)

% Choose default command line output for ODBCTool
handles.output = hObject;

set(get(handles.panel2, 'children'), 'enable', 'off');
handles.conn = [];
handles.db = [];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ODBCTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ODBCTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function btn_createds_Callback(hObject, eventdata, handles)
system('%SystemRoot%\system32\odbcad32.exe');

function lb_ds_Callback(hObject, eventdata, handles)


function lb_ds_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function btn_connODBC_Callback(hObject, eventdata, handles)
if get(hObject, 'value')%按下【ODBC连接】按钮
    data = get(handles.lb_ds, {'String', 'value'});
    if iscellstr(data{1}) %多个数据源
        handles.conn = database(data{1}{data{2}}, get(handles.username, 'string'), ...
            get(handles.password, 'string'));
    elseif ischar(data{1})
        handles.conn = database(data{1}, get(handles.username, 'string'), ...
            get(handles.password, 'string'));
    else
        return;
    end
    if isconnection(handles.conn)
        set(handles.tip, 'string', 'ODBC连接成功！');
        set(get(handles.panel2, 'children'), 'enable', 'on');            
        set(handles.listbox_db, 'string', catalogs(handles.conn));%query all db names            
    else
        set(handles.tip, 'string', 'ODBC连接失败！');
    end        
elseif ~isempty(handles.conn)
    handles.conn = [];
    set(handles.tip, 'string', '关闭连接成功。');
    set(get(handles.panel2, 'children'), 'enable', 'off');
end
guidata(hObject, handles);

function edit1_Callback(hObject, eventdata, handles)


function edit1_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit2_Callback(hObject, eventdata, handles)


function edit2_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)


function edit3_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function btn_sql_Callback(hObject, eventdata, handles)
if isempty(handles.conn)
    return;
end
sql1 = get(handles.sql, 'string');
a = sql1';
sql = a(:)';
curs = exec(handles.conn, sql);
curs2 = fetch(curs);
set(handles.datatable, 'data', curs2.Data,...
    'columnname', regexp(columnnames(curs2), '\w+', 'match'));
close(curs2);
close(curs);


function listbox_db_Callback(hObject, eventdata, handles)
if ~strcmp('open', get(gcf, 'selectionType'))
    data = get(hObject, {'string', 'value'});
    handles.db = data{1}{data{2}};
    set(handles.datatable, 'data', tables(handles.conn, handles.db),...
        'ColumnName', {'表名','表类型'});
    guidata(hObject, handles);
end

function listbox_db_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function btn_setdb_Callback(hObject, eventdata, handles)


function sql_Callback(hObject, eventdata, handles)


function sql_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function btn_queryds_Callback(hObject, eventdata, handles)
set(handles.lb_ds, 'string', getdatasources);


function username_Callback(hObject, eventdata, handles)


function username_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function password_Callback(hObject, eventdata, handles)


function password_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function password_KeyPressFcn(hObject, eventdata, handles)


function popmenu_db_Callback(hObject, eventdata, handles)


function popmenu_db_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function figure1_CloseRequestFcn(hObject, eventdata, handles)
if ~isempty(handles.conn)
    close(handles.conn);
end
% Hint: delete(hObject) closes the figure
delete(hObject);
