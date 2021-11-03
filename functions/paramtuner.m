function varargout = paramtuner(varargin)
% PARAMTUNER MATLAB code for paramtuner.fig
%      PARAMTUNER, by itself, creates a new PARAMTUNER or raises the existing
%      singleton*.
%
%      H = PARAMTUNER returns the handle to a new PARAMTUNER or the handle to
%      the existing singleton*.
%
%      PARAMTUNER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMTUNER.M with the given input arguments.
%
%      PARAMTUNER('Property','Value',...) creates a new PARAMTUNER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before paramtuner_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to paramtuner_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help paramtuner

% Last Modified by GUIDE v2.5 28-Sep-2021 16:28:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @paramtuner_OpeningFcn, ...
                   'gui_OutputFcn',  @paramtuner_OutputFcn, ...
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

% --- Executes just before paramtuner is made visible.
function paramtuner_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to paramtuner (see VARARGIN)

% Choose default command line output for paramtuner
handles.output = hObject;

myPath = fileparts(mfilename('fullpath'));
dataPath = fullfile(fileparts(myPath),'data');  
load(fullfile(dataPath,'time_series_all_channels.mat'))
handles.Time = Time;
handles.TimeDomainAligned = TimeDomainAligned;
load(fullfile(dataPath,'timeInformation.mat'),'timeInfo');
handles.infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;
% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using paramtuner.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

pushbutton2_Callback(hObject, eventdata, handles);



% UIWAIT makes paramtuner wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = paramtuner_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
fValues = linspace(0,250,1e3);
[F,FBL] = fittingmodel;
handles.params(10) = str2num(get(hObject,'String'));
axes(handles.axes1);
    cla; hold off;
    plot(handles.xData,handles.yData,'.-');
    hold on;
    hPlot1 = plot(fValues,F(handles.params,fValues));
    hPlot2 = plot(fValues,FBL(handles.params,fValues),'linewidth',1);
    xlim([min(handles.xData),max(handles.xData)]);
    set(gca,'xscale','log');
    gcaformat;
    xlabel('Frequency (Hz)');
guidata(hObject,handles);


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

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
fValues = linspace(0,250,1e3);
[F,FBL] = fittingmodel;
handles.params(11) = str2num(get(hObject,'String'));
axes(handles.axes1);
    cla; hold off;
    plot(handles.xData,handles.yData,'.-');
    hold on;
    hPlot1 = plot(fValues,F(handles.params,fValues));
    hPlot2 = plot(fValues,FBL(handles.params,fValues));
    xlim([min(handles.xData),max(handles.xData)]);
    set(gca,'xscale','log');
    gcaformat;
    xlabel('Frequency (Hz)');
guidata(hObject,handles);


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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
fValues = linspace(0,250,1e3);
[F,FBL] = fittingmodel;
handles.params(12) = str2num(get(hObject,'String'));
axes(handles.axes1);
    cla; hold off;
    plot(handles.xData,handles.yData,'.-');
    hold on;
    hPlot1 = plot(fValues,F(handles.params,fValues));
    hPlot2 = plot(fValues,FBL(handles.params,fValues));
    xlim([min(handles.xData),max(handles.xData)]);
    set(gca,'xscale','log');
    gcaformat;
    xlabel('Frequency (Hz)');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
fValues = linspace(0,250,1e3);
[F,FBL] = fittingmodel;
handles.params(13) = str2num(get(hObject,'String'));
axes(handles.axes1);
    cla; hold off;
    plot(handles.xData,handles.yData,'.-');
    hold on;
    hPlot1 = plot(fValues,F(handles.params,fValues));
    hPlot2 = plot(fValues,FBL(handles.params,fValues));
    xlim([min(handles.xData),max(handles.xData)]);
    set(gca,'xscale','log');
    gcaformat;
    xlabel('Frequency (Hz)');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
fValues = linspace(0,250,1e3);
[F,FBL] = fittingmodel;
handles.params(14) = str2num(get(hObject,'String'));
axes(handles.axes1);
    cla; hold off;
    plot(handles.xData,handles.yData,'.-');
    hold on;
    hPlot1 = plot(fValues,F(handles.params,fValues));
    hPlot2 = plot(fValues,FBL(handles.params,fValues));
    xlim([min(handles.xData),max(handles.xData)]);
    set(gca,'xscale','log');
    gcaformat;
    xlabel('Frequency (Hz)');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
fValues = linspace(0,250,1e3);
[F,FBL] = fittingmodel;
handles.params(15) = str2num(get(hObject,'String'));
axes(handles.axes1);
    cla; hold off;
    plot(handles.xData,handles.yData,'.-');
    hold on;
    hPlot1 = plot(fValues,F(handles.params,fValues));
    hPlot2 = plot(fValues,FBL(handles.params,fValues));
    xlim([min(handles.xData),max(handles.xData)]);
    set(gca,'xscale','log');
    gcaformat;
    xlabel('Frequency (Hz)');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(sum([strcmp(fieldnames(handles),'params')]))
    handles.params
   handles.params =  tunerGA(handles.popupmenu2.Value,~handles.checkbox1.Value,'N',50,'M',20,'fig',handles.axes1,'StartPoint',handles.params,'coolingOnset',1);
else
    handles.params =  tunerGA(handles.popupmenu2.Value,~handles.checkbox1.Value,'N',50,'M',20,'fig',handles.axes1);
end
handles.edit1.String = num2str(handles.params(10));
handles.edit2.String = num2str(handles.params(11));
handles.edit3.String = num2str(handles.params(12));
handles.edit4.String = num2str(handles.params(13));
handles.edit5.String = num2str(handles.params(14));
handles.edit6.String = num2str(handles.params(15));
guidata(hObject, handles);



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.params(1:9) =  fitGAsinglePeaks(handles.popupmenu2.Value,~handles.checkbox1.Value,'N',100,'M',60,'fig',handles.axes1,'StartPoint',handles.params);
guidata(hObject, handles);


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

pts = handles.popupmenu2.Value;
baseline = ~handles.checkbox1.Value;

[freq,time,psd] = eegfft(handles.Time,handles.TimeDomainAligned(:,2,pts),2,0.25);
if(baseline)
    it = interp1(time,1:length(time),handles.infusionTime(pts),'nearest');
    p = log(nanmedian(psd(:,1:it),2));
else
    it = interp1(time,1:length(time),0,'nearest');
    p = log(nanmedian(psd(:,it:end),2));
end
[xData, yData] = prepareCurveData( freq, p );
yData(xData>250) = [];
xData(xData>250) = [];
w = ones(size(xData));
for i = 1:floor(max(xData)/60)
    idcs = find(and(xData>60*i-5,xData<60*i+5));
    w(idcs) = 0;
end
idcs = find(w==1);
yData = interp1(xData(idcs),yData(idcs),xData,'linear');

handles.xData = xData;
handles.yData = yData;

fValues = linspace(0,250,1e3);
[F,FBL] = fittingmodel;
axes(handles.axes1);
    cla; hold off;
    plot(handles.xData,handles.yData,'.-');
    hold on;
    hPlot1 = plot(fValues,F(handles.params,fValues));
    hPlot2 = plot(fValues,FBL(handles.params,fValues));
    xlim([min(handles.xData),max(handles.xData)]);
    set(gca,'xscale','log');
    gcaformat;
    xlabel('Frequency (Hz)');
guidata(hObject,handles);
guidata(hObject, handles);

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

pts = handles.popupmenu2.Value;
baseline = ~handles.checkbox1.Value;

[freq,time,psd] = eegfft(handles.Time,handles.TimeDomainAligned(:,2,pts),2,0.25);
if(baseline)
    it = interp1(time,1:length(time),handles.infusionTime(pts),'nearest');
    p = log(nanmedian(psd(:,1:it),2));
else
    it = interp1(time,1:length(time),0,'nearest');
    p = log(nanmedian(psd(:,it:end),2));
end
[xData, yData] = prepareCurveData( freq, p );
yData(xData>250) = [];
xData(xData>250) = [];
w = ones(size(xData));
for i = 1:floor(max(xData)/60)
    idcs = find(and(xData>60*i-5,xData<60*i+5));
    w(idcs) = 0;
end
idcs = find(w==1);
yData = interp1(xData(idcs),yData(idcs),xData,'linear');

handles.xData = xData;
handles.yData = yData;

fValues = linspace(0,250,1e3);
[F,FBL] = fittingmodel;
axes(handles.axes1);
    cla; hold off;
    plot(handles.xData,handles.yData,'.-');
    hold on;
    hPlot1 = plot(fValues,F(handles.params,fValues));
    hPlot2 = plot(fValues,FBL(handles.params,fValues));
    xlim([min(handles.xData),max(handles.xData)]);
    set(gca,'xscale','log');
    gcaformat;
    xlabel('Frequency (Hz)');
guidata(hObject,handles);
guidata(hObject, handles);

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

round(handles.params(:),3,'significant')
