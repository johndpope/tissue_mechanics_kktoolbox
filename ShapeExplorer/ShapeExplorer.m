function varargout = ShapeExplorer(varargin)
% SHAPEEXPLORER M-file for ShapeExplorer.fig
%      SHAPEEXPLORER, by itself, creates a new SHAPEEXPLORER or raises the existing
%      singleton*.
%
%      H = SHAPEEXPLORER returns the handle to a new SHAPEEXPLORER or the handle to
%      the existing singleton*.
%
%      SHAPEEXPLORER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHAPEEXPLORER.M with the given input arguments.
%
%      SHAPEEXPLORER('Property','Value',...) creates a new SHAPEEXPLORER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ShapeExplorer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ShapeExplorer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ShapeExplorer

% Last Modified by GUIDE v2.5 21-Mar-2012 10:55:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ShapeExplorer_OpeningFcn, ...
                   'gui_OutputFcn',  @ShapeExplorer_OutputFcn, ...
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


% --- Executes just before ShapeExplorer is made visible.
function ShapeExplorer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ShapeExplorer (see VARARGIN)

%%%% load some preliminary data
disp('initializing');
%%% define the basis
nico = 3;
gdim = 60;
L_max = 12;
X_o = zeros(3*(L_max + 1)^2);

[X,C]=surface_mesh.sphere_mesh_gen(nico);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
[x y z] = kk_sph2cart(t,p,1);
%%% generate the new basis
[L, K] = shp_surface.indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
Y_LK  = zeros(N, M, 'single');
for S = 1:length(L),
    Y_LK(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;
%%%%%%
%%% for transition to objects generate a basis object and s
basis = sh_basis(L_max, gdim);
s = shp_surface(basis);
s.X_o = X_o;
s = update(s);
handles.ud.s = s;
handles.ud.Y_LK = Y_LK;
handles.ud.C = C;
handles.ud.X_o = X_o;
handles.ud.fac = 1;     % factor by which to multiply the clks to obtain real units
handles.ud.coordinate = 1;
handles.ud.gdim = gdim;

% Choose default command line output for sc
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
initialize_Callback(hObject,eventdata,handles);
disp('done');

% UIWAIT makes ShapeExplorer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ShapeExplorer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 1;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
axis off; 

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 2;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in initialize.
function initialize_Callback(hObject, eventdata, handles)
% hObject    handle to initialize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ud.X_o = zeros(size(handles.ud.X_o));
[xclks yclks zclks] = shp_surface.get_xyz_clks(handles.ud.X_o);
xclks(4) = -1;yclks(2) = -1;zclks(3) = 0.3;
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];
guidata(hObject, handles);% Update handles structure
shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;
synchronize_slider_ranges(handles);
synchronize(handles);
set(handles.figure1,'Name','untitled');



% --- Executes on selection change in choose_coordinate.
function choose_coordinate_Callback(hObject, eventdata, handles)
% hObject    handle to choose_coordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns choose_coordinate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from choose_coordinate
handles.ud.coordinate = get(hObject,'Value');
synchronize_slider_ranges(handles);
synchronize(handles);
guidata(hObject, handles);% Update handles structure

% --- Executes during object creation, after setting all properties.
function choose_coordinate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to choose_coordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 3;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 4;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 5;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 6;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 7;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 8;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 9;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider10_Callback(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 10;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_slider_ranges(handles)
% set the values of all slider min max values to match the coefficients of the current
% coordinate
for ix  = 1:169%length(coord_clks),
    str = sprintf(' set(handles.slider%d, ''Min'',%.2f);',ix,-1);eval(str);
    str = sprintf(' set(handles.slider%d, ''Max'',%.2f);',ix,1);eval(str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function synchronize_slider_ranges(handles)
% set the values of all slider min max values to match the coefficients of the current
% coordinate
initialize_slider_ranges(handles);
clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if     handles.ud.coordinate ==1,coord_clks = xclks;
elseif handles.ud.coordinate ==2,coord_clks = yclks;
elseif handles.ud.coordinate ==3,coord_clks = zclks;end
for ix  = 1:169%length(coord_clks),
    str = sprintf('if (get(handles.slider%d, ''Min'')> %.2f), set(handles.slider%d, ''Min'',%.2f);end;',ix, coord_clks(ix), ix,coord_clks(ix));eval(str);
    str = sprintf('if (get(handles.slider%d, ''Max'')< %.2f), set(handles.slider%d, ''Max'',%.2f);end;',ix, coord_clks(ix), ix,coord_clks(ix));eval(str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function synchronize(handles)
% set the values of all sliders to match the coefficients of the current
% coordinate
clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if     handles.ud.coordinate ==1,coord_clks = xclks;
elseif handles.ud.coordinate ==2,coord_clks = yclks;
elseif handles.ud.coordinate ==3,coord_clks = zclks;end
for ix  = 1:169%length(coord_clks),
    str = sprintf('set(handles.slider%d, ''Value'',%.2f);',ix,coord_clks(ix));eval(str);
end



% --- Executes on slider movement.
function slider11_Callback(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 11;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider12_Callback(hObject, eventdata, handles)
% hObject    handle to slider12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 12;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider13_Callback(hObject, eventdata, handles)
% hObject    handle to slider13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 13;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider14_Callback(hObject, eventdata, handles)
% hObject    handle to slider14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 14;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider15_Callback(hObject, eventdata, handles)
% hObject    handle to slider15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 15;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider16_Callback(hObject, eventdata, handles)
% hObject    handle to slider16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 16;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider17_Callback(hObject, eventdata, handles)
% hObject    handle to slider17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 17;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider18_Callback(hObject, eventdata, handles)
% hObject    handle to slider18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 18;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider20_Callback(hObject, eventdata, handles)
% hObject    handle to slider20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 20;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider21_Callback(hObject, eventdata, handles)
% hObject    handle to slider21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 21;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider22_Callback(hObject, eventdata, handles)
% hObject    handle to slider22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 22;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider23_Callback(hObject, eventdata, handles)
% hObject    handle to slider23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 23;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider19_Callback(hObject, eventdata, handles)
% hObject    handle to slider19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 19;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on slider movement.
function slider24_Callback(hObject, eventdata, handles)
% hObject    handle to slider24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 24;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider25_Callback(hObject, eventdata, handles)
% hObject    handle to slider25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 25;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider26_Callback(hObject, eventdata, handles)
% hObject    handle to slider26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 26;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider27_Callback(hObject, eventdata, handles)
% hObject    handle to slider27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 27;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider28_Callback(hObject, eventdata, handles)
% hObject    handle to slider28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 28;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider29_Callback(hObject, eventdata, handles)
% hObject    handle to slider29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 29;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider30_Callback(hObject, eventdata, handles)
% hObject    handle to slider30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 30;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider31_Callback(hObject, eventdata, handles)
% hObject    handle to slider31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 31;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider32_Callback(hObject, eventdata, handles)
% hObject    handle to slider32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 32;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider34_Callback(hObject, eventdata, handles)
% hObject    handle to slider25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 34;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider33_Callback(hObject, eventdata, handles)
% hObject    handle to slider33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 33;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider35_Callback(hObject, eventdata, handles)
% hObject    handle to slider35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 35;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider36_Callback(hObject, eventdata, handles)
% hObject    handle to slider36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 36;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in import.
function import_Callback(hObject, eventdata, handles)
% hObject    handle to import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterindex] = uigetfile('*.mat', 'Select a shape MAT-file to explore');
if isequal(filename,0) || isequal(pathname,0)
%     disp('User pressed cancel')
else
    load(fullfile(pathname, filename));
    set(handles.figure1,'Name',fullfile(pathname, filename));
    L_max = sqrt(length(handles.ud.X_o)/3)-1;
    %%% crop the shape vector according to precalculated basis set
    if exist('X_o','var'),
        clks = X_o;
        nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
    elseif exist('xclks','var'),
        clks = [xclks(:)' yclks(:)' zclks(:)'];
    else
        clks = handles.ud.X_o;
        nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
    end
    handles.ud.fac = max(abs(clks));
    disp(handles.ud.fac);
%     xclks = xclks/max(abs(clks));yclks = yclks/max(abs(clks));zclks = zclks/(max(abs(clks)));
    if length(xclks)>(L_max+1)^2, xclks = xclks(1:(L_max+1)^2);yclks = yclks(1:(L_max+1)^2);zclks = zclks(1:(L_max+1)^2);end
    if length(xclks)<(L_max+1)^2, xclks((L_max+1)^2) = 0;yclks((L_max+1)^2) = 0;zclks((L_max+1)^2) = 0;end
    handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];
    guidata(hObject, handles);
    shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;axis tight;axis manual;
    synchronize_slider_ranges(handles);
    synchronize(handles);
end


% --- Executes on button press in center.
function center_Callback(hObject, eventdata, handles)
% hObject    handle to center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axis tight;axis manual;


% --- Executes on slider movement.
function slider37_Callback(hObject, eventdata, handles)
% hObject    handle to slider37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 37;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider38_Callback(hObject, eventdata, handles)
% hObject    handle to slider38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 38;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider39_Callback(hObject, eventdata, handles)
% hObject    handle to slider39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 39;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider40_Callback(hObject, eventdata, handles)
% hObject    handle to slider40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 40;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider42_Callback(hObject, eventdata, handles)
% hObject    handle to slider42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 42;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider43_Callback(hObject, eventdata, handles)
% hObject    handle to slider43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 43;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider41_Callback(hObject, eventdata, handles)
% hObject    handle to slider41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 41;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider44_Callback(hObject, eventdata, handles)
% hObject    handle to slider44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 44;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider45_Callback(hObject, eventdata, handles)
% hObject    handle to slider45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 45;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider46_Callback(hObject, eventdata, handles)
% hObject    handle to slider46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 46;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider47_Callback(hObject, eventdata, handles)
% hObject    handle to slider47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 47;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider48_Callback(hObject, eventdata, handles)
% hObject    handle to slider48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 48;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider49_Callback(hObject, eventdata, handles)
% hObject    handle to slider49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

indx = 49;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on slider movement.
function slider51_Callback(hObject, eventdata, handles)
% hObject    handle to slider51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 51;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider52_Callback(hObject, eventdata, handles)
% hObject    handle to slider52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 52;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider53_Callback(hObject, eventdata, handles)
% hObject    handle to slider53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 53;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider53_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider55_Callback(hObject, eventdata, handles)
% hObject    handle to slider55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 55;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider56_Callback(hObject, eventdata, handles)
% hObject    handle to slider56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 56;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider56_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider54_Callback(hObject, eventdata, handles)
% hObject    handle to slider54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 54;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider54_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider57_Callback(hObject, eventdata, handles)
% hObject    handle to slider57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 57;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider57_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider58_Callback(hObject, eventdata, handles)
% hObject    handle to slider58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 58;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider58_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider59_Callback(hObject, eventdata, handles)
% hObject    handle to slider59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 59;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider59_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider73_Callback(hObject, eventdata, handles)
% hObject    handle to slider73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 73;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider73_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider61_Callback(hObject, eventdata, handles)
% hObject    handle to slider61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 61;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider61_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider62_Callback(hObject, eventdata, handles)
% hObject    handle to slider62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 62;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider62_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider50_Callback(hObject, eventdata, handles)
% hObject    handle to slider50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 50;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure

% --- Executes during object creation, after setting all properties.
function slider50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider60_Callback(hObject, eventdata, handles)
% hObject    handle to slider60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 60;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider60_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider63_Callback(hObject, eventdata, handles)
% hObject    handle to slider63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 63;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider63_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider64_Callback(hObject, eventdata, handles)
% hObject    handle to slider64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 64;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider64_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X_o = handles.ud.X_o;
uisave('X_o');



% --- Executes on slider movement.
function slider66_Callback(hObject, eventdata, handles)
% hObject    handle to slider66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 66;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider67_Callback(hObject, eventdata, handles)
% hObject    handle to slider67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 67;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider67_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider68_Callback(hObject, eventdata, handles)
% hObject    handle to slider68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 68;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider68_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider70_Callback(hObject, eventdata, handles)
% hObject    handle to slider70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 70;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider70_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider71_Callback(hObject, eventdata, handles)
% hObject    handle to slider71 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 71;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider71_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider71 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider69_Callback(hObject, eventdata, handles)
% hObject    handle to slider69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 69;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider69_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider72_Callback(hObject, eventdata, handles)
% hObject    handle to slider72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 72;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider72_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider85_Callback(hObject, eventdata, handles)
% hObject    handle to slider73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 85;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider85_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider74_Callback(hObject, eventdata, handles)
% hObject    handle to slider74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 74;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider74_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider75_Callback(hObject, eventdata, handles)
% hObject    handle to slider75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 75;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider75_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider76_Callback(hObject, eventdata, handles)
% hObject    handle to slider76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 76;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider76_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider77_Callback(hObject, eventdata, handles)
% hObject    handle to slider77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 77;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider77_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider65_Callback(hObject, eventdata, handles)
% hObject    handle to slider65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 65;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider65_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider78_Callback(hObject, eventdata, handles)
% hObject    handle to slider78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

indx = 78;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider78_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider79_Callback(hObject, eventdata, handles)
% hObject    handle to slider79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 79;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider79_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider80_Callback(hObject, eventdata, handles)
% hObject    handle to slider80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 80;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider80_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider81_Callback(hObject, eventdata, handles)
% hObject    handle to slider81 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 81;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider81_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider81 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider83_Callback(hObject, eventdata, handles)
% hObject    handle to slider83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

indx = 83;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider83_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider84_Callback(hObject, eventdata, handles)
% hObject    handle to slider84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 84;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider84_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider97_Callback(hObject, eventdata, handles)
% hObject    handle to slider85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 97;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider97_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider87_Callback(hObject, eventdata, handles)
% hObject    handle to slider87 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 87;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider87_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider87 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider99_Callback(hObject, eventdata, handles)
% hObject    handle to slider99 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 99;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider99_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider99 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider86_Callback(hObject, eventdata, handles)
% hObject    handle to slider86 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 86;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider86_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider86 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider89_Callback(hObject, eventdata, handles)
% hObject    handle to slider89 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 89;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider89_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider89 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider90_Callback(hObject, eventdata, handles)
% hObject    handle to slider90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 90;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider90_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider91_Callback(hObject, eventdata, handles)
% hObject    handle to slider91 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 91;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider91_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider91 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider92_Callback(hObject, eventdata, handles)
% hObject    handle to slider92 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 92;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider92_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider92 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider93_Callback(hObject, eventdata, handles)
% hObject    handle to slider93 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 93;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider93_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider93 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider94_Callback(hObject, eventdata, handles)
% hObject    handle to slider94 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 94;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider94_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider94 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider82_Callback(hObject, eventdata, handles)
% hObject    handle to slider82 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 82;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider82_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider82 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider95_Callback(hObject, eventdata, handles)
% hObject    handle to slider95 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 95;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider95_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider95 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider96_Callback(hObject, eventdata, handles)
% hObject    handle to slider96 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 96;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider96_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider96 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on slider movement.
function slider98_Callback(hObject, eventdata, handles)
% hObject    handle to slider98 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 98;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider98_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider98 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider100_Callback(hObject, eventdata, handles)
% hObject    handle to slider100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 100;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider100_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on slider movement.
function slider88_Callback(hObject, eventdata, handles)
% hObject    handle to slider88 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 88;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure



% --- Executes during object creation, after setting all properties.
function slider88_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider88 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on slider movement.
function slider114_Callback(hObject, eventdata, handles)
% hObject    handle to slider114 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 114;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider114_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider114 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider115_Callback(hObject, eventdata, handles)
% hObject    handle to slider115 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 115;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider115_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider115 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider116_Callback(hObject, eventdata, handles)
% hObject    handle to slider116 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 116;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider116_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider116 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider117_Callback(hObject, eventdata, handles)
% hObject    handle to slider117 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 117;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider117_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider117 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider118_Callback(hObject, eventdata, handles)
% hObject    handle to slider118 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 118;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider118_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider118 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider119_Callback(hObject, eventdata, handles)
% hObject    handle to slider119 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 119;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider119_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider119 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider120_Callback(hObject, eventdata, handles)
% hObject    handle to slider120 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 120;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider120_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider120 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider121_Callback(hObject, eventdata, handles)
% hObject    handle to slider121 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 121;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider121_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider121 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider122_Callback(hObject, eventdata, handles)
% hObject    handle to slider122 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 122;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider122_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider122 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider123_Callback(hObject, eventdata, handles)
% hObject    handle to slider123 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 123;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider123_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider123 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider124_Callback(hObject, eventdata, handles)
% hObject    handle to slider124 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 124;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider124_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider124 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider125_Callback(hObject, eventdata, handles)
% hObject    handle to slider125 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 125;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider125_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider125 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider126_Callback(hObject, eventdata, handles)
% hObject    handle to slider126 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 126;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider126_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider126 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider127_Callback(hObject, eventdata, handles)
% hObject    handle to slider127 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 127;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider127_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider127 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider128_Callback(hObject, eventdata, handles)
% hObject    handle to slider128 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 128;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider128_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider128 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider129_Callback(hObject, eventdata, handles)
% hObject    handle to slider129 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 129;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider129_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider129 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider130_Callback(hObject, eventdata, handles)
% hObject    handle to slider130 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 130;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider130_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider130 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider131_Callback(hObject, eventdata, handles)
% hObject    handle to slider131 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 131;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider131_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider131 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider132_Callback(hObject, eventdata, handles)
% hObject    handle to slider132 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 132;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider132_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider132 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider133_Callback(hObject, eventdata, handles)
% hObject    handle to slider133 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 133;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider133_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider133 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider134_Callback(hObject, eventdata, handles)
% hObject    handle to slider134 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 134;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider134_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider134 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider135_Callback(hObject, eventdata, handles)
% hObject    handle to slider135 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 135;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider135_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider135 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider136_Callback(hObject, eventdata, handles)
% hObject    handle to slider136 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 136;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider136_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider136 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider137_Callback(hObject, eventdata, handles)
% hObject    handle to slider137 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 137;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider137_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider137 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider138_Callback(hObject, eventdata, handles)
% hObject    handle to slider138 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 138;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider138_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider138 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider139_Callback(hObject, eventdata, handles)
% hObject    handle to slider139 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 139;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider139_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider139 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider140_Callback(hObject, eventdata, handles)
% hObject    handle to slider140 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 140;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider140_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider140 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider141_Callback(hObject, eventdata, handles)
% hObject    handle to slider141 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 141;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider141_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider141 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider142_Callback(hObject, eventdata, handles)
% hObject    handle to slider142 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 142;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider142_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider142 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider143_Callback(hObject, eventdata, handles)
% hObject    handle to slider143 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 143;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider143_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider143 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider144_Callback(hObject, eventdata, handles)
% hObject    handle to slider144 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 144;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider144_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider144 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider145_Callback(hObject, eventdata, handles)
% hObject    handle to slider145 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 145;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider145_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider145 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider146_Callback(hObject, eventdata, handles)
% hObject    handle to slider146 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 146;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider146_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider146 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider147_Callback(hObject, eventdata, handles)
% hObject    handle to slider147 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 147;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider147_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider147 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider148_Callback(hObject, eventdata, handles)
% hObject    handle to slider148 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 148;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider148_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider148 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider149_Callback(hObject, eventdata, handles)
% hObject    handle to slider149 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 149;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider149_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider149 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider150_Callback(hObject, eventdata, handles)
% hObject    handle to slider150 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 150;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider150_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider150 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider151_Callback(hObject, eventdata, handles)
% hObject    handle to slider151 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 151;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider151_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider151 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider152_Callback(hObject, eventdata, handles)
% hObject    handle to slider152 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 152;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider152_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider152 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider153_Callback(hObject, eventdata, handles)
% hObject    handle to slider153 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 153;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider153_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider153 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider154_Callback(hObject, eventdata, handles)
% hObject    handle to slider154 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 154;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider154_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider154 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider155_Callback(hObject, eventdata, handles)
% hObject    handle to slider155 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 155;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider155_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider155 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider101_Callback(hObject, eventdata, handles)
% hObject    handle to slider101 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 101;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider101_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider101 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider102_Callback(hObject, eventdata, handles)
% hObject    handle to slider102 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 102;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider102_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider102 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider156_Callback(hObject, eventdata, handles)
% hObject    handle to slider156 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 156;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider156_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider156 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider157_Callback(hObject, eventdata, handles)
% hObject    handle to slider157 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 157;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider157_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider157 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider158_Callback(hObject, eventdata, handles)
% hObject    handle to slider158 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 158;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider158_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider158 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider159_Callback(hObject, eventdata, handles)
% hObject    handle to slider159 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 159;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider159_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider159 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider160_Callback(hObject, eventdata, handles)
% hObject    handle to slider160 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 160;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider160_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider160 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider161_Callback(hObject, eventdata, handles)
% hObject    handle to slider161 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 161;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider161_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider161 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider162_Callback(hObject, eventdata, handles)
% hObject    handle to slider162 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 162;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider162_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider162 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider163_Callback(hObject, eventdata, handles)
% hObject    handle to slider163 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 163;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider163_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider163 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider164_Callback(hObject, eventdata, handles)
% hObject    handle to slider164 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 164;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider164_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider164 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider165_Callback(hObject, eventdata, handles)
% hObject    handle to slider165 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 165;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider165_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider165 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider166_Callback(hObject, eventdata, handles)
% hObject    handle to slider166 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 166;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider166_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider166 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider167_Callback(hObject, eventdata, handles)
% hObject    handle to slider167 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 167;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider167_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider167 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider168_Callback(hObject, eventdata, handles)
% hObject    handle to slider168 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 168;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider168_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider168 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider169_Callback(hObject, eventdata, handles)
% hObject    handle to slider169 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 169;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider169_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider169 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider103_Callback(hObject, eventdata, handles)
% hObject    handle to slider103 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 103;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider103_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider103 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider104_Callback(hObject, eventdata, handles)
% hObject    handle to slider104 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 104;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider104_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider104 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider105_Callback(hObject, eventdata, handles)
% hObject    handle to slider105 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 105;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider105_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider105 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider106_Callback(hObject, eventdata, handles)
% hObject    handle to slider106 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 106;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider106_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider106 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider107_Callback(hObject, eventdata, handles)
% hObject    handle to slider107 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 107;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider107_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider107 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider108_Callback(hObject, eventdata, handles)
% hObject    handle to slider108 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 108;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider108_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider108 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider109_Callback(hObject, eventdata, handles)
% hObject    handle to slider109 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 109;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider109_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider109 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider110_Callback(hObject, eventdata, handles)
% hObject    handle to slider110 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 110;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider110_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider110 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider111_Callback(hObject, eventdata, handles)
% hObject    handle to slider111 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 111;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider111_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider111 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider112_Callback(hObject, eventdata, handles)
% hObject    handle to slider112 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 112;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider112_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider112 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider113_Callback(hObject, eventdata, handles)
% hObject    handle to slider113 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = 113;clks  = handles.ud.X_o;nc = length(clks)/3;xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks = clks(2*nc+1:end);
if handles.ud.coordinate ==1, xclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==2, yclks(indx) = get(hObject, 'Value');elseif handles.ud.coordinate ==3,zclks(indx) = get(hObject, 'Value');end
handles.ud.X_o = [xclks(:)' yclks(:)' zclks(:)'];shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure


% --- Executes during object creation, after setting all properties.
function slider113_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider113 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in Calc_props.
function Calc_props_Callback(hObject, eventdata, handles)
% hObject    handle to Calc_props (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [Area,V,v_red,t,p,X,C,Y_LK]=shape_explorer_plot_sh(handles.ud.X_o * handles.ud.fac,handles.ud.gdim, handles.axes1);axis tight;axis manual;
[xclks yclks zclks] = shp_surface.get_xyz_clks(handles.ud.X_o);
X = handles.ud.Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
[A, V, v, F_areas, h, H, Eb, da] = shape_explorer_triangulated_props(X, handles.ud.C, 1);
%% test for self-intersection
handles.ud.s.X_o = handles.ud.X_o;
s = update(handles.ud.s);
si = self_intersection(s);
disp(['Self-intersection: ' num2str(si)]);
rotate3d on;




function Lmax_Callback(hObject, eventdata, handles)
% hObject    handle to Lmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lmax as text
%        str2double(get(hObject,'String')) returns contents of Lmax as a double
newLmax = str2double(get(hObject,'String'));

% clks = handles.ud.X_o;
% nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
% handles.ud.Lmax = newLmax;

% --- Executes during object creation, after setting all properties.
function Lmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gdim_Callback(hObject, eventdata, handles)
% hObject    handle to gdim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gdim as text
%        str2double(get(hObject,'String')) returns contents of gdim as a double
handles.ud.gdim = str2num(get(hObject,'String'));
shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);axis tight;axis manual;guidata(hObject, handles);% Update handles structure

% --- Executes during object creation, after setting all properties.
function gdim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gdim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in rotation_invariant.
function rotation_invariant_Callback(hObject, eventdata, handles)
% hObject    handle to rotation_invariant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ud.s.X_o = handles.ud.X_o;
handles.ud.s = update(handles.ud.s);
disp('Calculating rotational invariant form ...');
handles.ud.s = r_inv(handles.ud.s);
disp('Done!');
handles.ud.X_o = handles.ud.s.X_o;
shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, handles.axes1);
axis tight;axis manual; drawnow;rotate3d on
    synchronize_slider_ranges(handles);
    synchronize(handles);
guidata(hObject, handles);% Update handles structure




% --- Executes on button press in sep_fig.
function sep_fig_Callback(hObject, eventdata, handles)
% hObject    handle to sep_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = figure;shape_explorer_plot_sh(handles.ud.X_o,handles.ud.Y_LK, handles.ud.C, gca(h));