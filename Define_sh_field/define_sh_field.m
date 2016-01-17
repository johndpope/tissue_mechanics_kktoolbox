function varargout = define_sh_field(varargin)
% DEFINE_SH_FIELD MATLAB code for define_sh_field.fig
%      DEFINE_SH_FIELD, by itself, creates a new DEFINE_SH_FIELD or raises the existing
%      singleton*.
%
%      H = DEFINE_SH_FIELD returns the handle to a new DEFINE_SH_FIELD or the handle to
%      the existing singleton*.
%
%      DEFINE_SH_FIELD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEFINE_SH_FIELD.M with the given input arguments.
%
%      DEFINE_SH_FIELD('Property','Value',...) creates a new DEFINE_SH_FIELD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before define_sh_field_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to define_sh_field_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help define_sh_field

% Last Modified by GUIDE v2.5 08-May-2012 12:12:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @define_sh_field_OpeningFcn, ...
                   'gui_OutputFcn',  @define_sh_field_OutputFcn, ...
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


% --- Executes just before define_sh_field is made visible.
function define_sh_field_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to define_sh_field (see VARARGIN)

   % make this the current axis
nico = 5;
d = dros_embryo('X_o_dros_embryo_04.mat');
s = sh_surface(d.L_max, d.basis);
[Y_LK C] = plotting_basis_gen(d,nico);
% Choose default command line output for define_sh_field
%% define the sh_surface object and assign it as output

handles.output = hObject;

% Update handles structure
handles.name = 'untitled';
handles.d = d;
handles.s = s;
handles.f = ones((d.L_max+1)^2,1);    % multiplicative factor for controling range and sensitivity for each clk
handles.Y_LK = Y_LK;
handles.C = C;
%% plot the drosophila embryo outline and enable 3D rotation
plot_dros(handles);

% Choose default command line output for define_sh_field
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes define_sh_field wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = define_sh_field_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    %% add the current field to the drosophila object
    handles.d = add_sf(handles.d,handles.name,handles.s);
    % Get default command line output from handles structure
    varargout{1} = handles.d;
    varargout{2} = handles.s;
    % The GUI is no longer waiting, just close it
    delete(hObject);
end
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% utility functions %%%%%%%%%%%%%%%%%%%
function plot_dros(handles)
obj = handles.d;

%% generate the shape and scalar field
X = handles.Y_LK(:,1:length(handles.d.xc))* [handles.d.xc(:) handles.d.yc(:) handles.d.zc(:)];
sf = handles.Y_LK(:,1:length(handles.s.xc))*handles.s.xc(:);
%% do the actual plotting
axes(handles.axes1);cla;patch('Vertices', X, 'Faces', handles.C,'FaceVertexCData',sf,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);daspect([1 1 1]);axis off;lighting gouraud;view(0,0);camlight;rotate3d
axes(handles.axes2);cla;patch('Vertices', X, 'Faces', handles.C,'FaceVertexCData',sf,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);daspect([1 1 1]);axis tight;view(0,-90);
axes(handles.axes3);cla;patch('Vertices', X, 'Faces', handles.C,'FaceVertexCData',sf,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);daspect([1 1 1]);axis tight;view(0,-270);

axes(handles.axes1);rotate3d;

function [Y_LK C] = plotting_basis_gen(obj,nico)
%% generate subdivisions of icosahedron
if nico > 6, nico = 6; disp(['Icosahedron subdivision: ' num2str(nico)]);end;
[X,C]=surface_mesh.sphere_mesh_gen(nico);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
% [x y z] = kk_sph2cart(t,p,1);
%% % generate the new basis
[L, K] = obj.indices_gen(1:(obj.L_max + 1)^2); M = length(L);N = length(t);
Y_LK  = zeros(N, M, 'double');
for S = 1:length(L),
    Y_LK(:,S) = obj.basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;
% --- Executes on button press in pushbutton_initialize.
function pushbutton_initialize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_initialize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.s.xc = zeros(size(handles.s.xc));
plot_dros(handles)

%%% set all sliders to zero
set(handles.slider00,'Value',0);
set(handles.slider1m1,'Value',0);
set(handles.slider10,'Value',0);
set(handles.slider11,'Value',0);
set(handles.slider2m2,'Value',0);
set(handles.slider2m1,'Value',0);
set(handles.slider20,'Value',0);
set(handles.slider21,'Value',0);
set(handles.slider22,'Value',0);
set(handles.slider3m3,'Value',0);
set(handles.slider3m2,'Value',0);
set(handles.slider3m1,'Value',0);
set(handles.slider30,'Value',0);
set(handles.slider31,'Value',0);
set(handles.slider32,'Value',0);
set(handles.slider33,'Value',0);
set(handles.slider4m4,'Value',0);
set(handles.slider4m3,'Value',0);
set(handles.slider4m2,'Value',0);
set(handles.slider4m1,'Value',0);
set(handles.slider40,'Value',0);
set(handles.slider41,'Value',0);
set(handles.slider42,'Value',0);
set(handles.slider43,'Value',0);
set(handles.slider44,'Value',0);
set(handles.slider5m5,'Value',0);
set(handles.slider5m4,'Value',0);
set(handles.slider5m3,'Value',0);
set(handles.slider5m2,'Value',0);
set(handles.slider5m1,'Value',0);
set(handles.slider50,'Value',0);
set(handles.slider51,'Value',0);
set(handles.slider52,'Value',0);
set(handles.slider53,'Value',0);
set(handles.slider54,'Value',0);
set(handles.slider55,'Value',0);
set(handles.slider6m6,'Value',0);
set(handles.slider6m5,'Value',0);
set(handles.slider6m4,'Value',0);
set(handles.slider6m3,'Value',0);
set(handles.slider6m2,'Value',0);
set(handles.slider6m1,'Value',0);
set(handles.slider60,'Value',0);
set(handles.slider61,'Value',0);
set(handles.slider62,'Value',0);
set(handles.slider63,'Value',0);
set(handles.slider64,'Value',0);
set(handles.slider65,'Value',0);
set(handles.slider66,'Value',0);

guidata(hObject,handles);



% --- Executes on slider movement.
function slider00_Callback(hObject, eventdata, handles)
% hObject    handle to slider00 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(211) = handles.f(2)*p; % 157 for L = 12, K = 0; 211 for L = 14, K = 0
plot_dros(handles);
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function slider00_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider00 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider1m1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(2) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider1m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1m1 (see GCBO)
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
p = get(hObject,'Value');
handles.s.xc(3) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider11_Callback(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(4) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
function slider2m2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(5) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider2m2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2m1_Callback(hObject, eventdata, handles)
% hObject    handle to slider2m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(6) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider2m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2m1 (see GCBO)
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
p = get(hObject,'Value');
handles.s.xc(7) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
p = get(hObject,'Value');
handles.s.xc(8) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
p = get(hObject,'Value');
handles.s.xc(9) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
function slider3m3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3m3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(10) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider3m3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3m3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3m2_Callback(hObject, eventdata, handles)
% hObject    handle to slider3m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(11) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider3m2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3m1_Callback(hObject, eventdata, handles)
% hObject    handle to slider3m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(12) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider3m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3m1 (see GCBO)
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
p = get(hObject,'Value');
handles.s.xc(13) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
p = get(hObject,'Value');
handles.s.xc(14) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
p = get(hObject,'Value');
handles.s.xc(15) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
function slider33_Callback(hObject, eventdata, handles)
% hObject    handle to slider33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(16) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
function slider4m4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4m4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(17) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider4m4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4m4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4m3_Callback(hObject, eventdata, handles)
% hObject    handle to slider4m3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(18) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider4m3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4m3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4m2_Callback(hObject, eventdata, handles)
% hObject    handle to slider4m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(19) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider4m2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4m1_Callback(hObject, eventdata, handles)
% hObject    handle to slider4m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(20) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider4m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4m1 (see GCBO)
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
p = get(hObject,'Value');
handles.s.xc(21) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
function slider41_Callback(hObject, eventdata, handles)
% hObject    handle to slider41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(22) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
function slider42_Callback(hObject, eventdata, handles)
% hObject    handle to slider42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');
handles.s.xc(23) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
p = get(hObject,'Value');
handles.s.xc(24) = handles.f(2)*p;
plot_dros(handles);
guidata(hObject,handles);

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
function slider44_Callback(hObject, eventdata, handles)
% hObject    handle to slider44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(25) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);

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
function slider5m5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5m5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(26) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider5m5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5m5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5m4_Callback(hObject, eventdata, handles)
% hObject    handle to slider5m4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(27) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider5m4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5m4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5m3_Callback(hObject, eventdata, handles)
% hObject    handle to slider5m3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(28) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider5m3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5m3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5m2_Callback(hObject, eventdata, handles)
% hObject    handle to slider5m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(29) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider5m2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5m1_Callback(hObject, eventdata, handles)
% hObject    handle to slider5m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(30) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider5m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5m1 (see GCBO)
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
p = get(hObject,'Value');handles.s.xc(31) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


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
function slider51_Callback(hObject, eventdata, handles)
% hObject    handle to slider51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(32) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


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
p = get(hObject,'Value');handles.s.xc(33) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


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
p = get(hObject,'Value');handles.s.xc(34) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


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
function slider54_Callback(hObject, eventdata, handles)
% hObject    handle to slider54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(35) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


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
function slider55_Callback(hObject, eventdata, handles)
% hObject    handle to slider55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(36) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


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
function slider6m6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6m6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(37) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider6m6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6m6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6m5_Callback(hObject, eventdata, handles)
% hObject    handle to slider6m5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(38) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider6m5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6m5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6m4_Callback(hObject, eventdata, handles)
% hObject    handle to slider6m4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(39) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider6m4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6m4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6m3_Callback(hObject, eventdata, handles)
% hObject    handle to slider6m3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(40) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider6m3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6m3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6m2_Callback(hObject, eventdata, handles)
% hObject    handle to slider6m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(41) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider6m2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6m1_Callback(hObject, eventdata, handles)
% hObject    handle to slider6m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(42) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider6m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6m1 (see GCBO)
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
p = get(hObject,'Value');handles.s.xc(43) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


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
function slider61_Callback(hObject, eventdata, handles)
% hObject    handle to slider61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(44) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


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
p = get(hObject,'Value');handles.s.xc(45) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


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
function slider63_Callback(hObject, eventdata, handles)
% hObject    handle to slider63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(46) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


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
p = get(hObject,'Value');handles.s.xc(47) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider64_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider64 (see GCBO)
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

p = get(hObject,'Value');handles.s.xc(48) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);

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
function slider66_Callback(hObject, eventdata, handles)
% hObject    handle to slider66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p = get(hObject,'Value');handles.s.xc(49) = handles.f(2)*p;plot_dros(handles);guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xc = handles.s.xc;
uisave({'xc'},'SF_field');

% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiopen('load');
%%% resize xc if necessary to be compatible with the current basis
lmax_in = sqrt(length(xc))-1;
trunc = (handles.d.L_max+1)^2;
if lmax_in<handles.d.L_max,
    xc(trunc) = 0;  % pads xc vector to the required length
else
    % else we need to pad handles.d (i.e. generate a new basis)
    warning('Should generate new basis to match field --- feature not programmed yet-- truncating field');
    xc = xc(1:trunc);
end
%%%%%%%%%%%%%%%%%%
handles.s.xc = xc;
s = handles.s;
plot_dros(handles)

%%% set all sliders to zero
set(handles.slider00,'Value',s.xc(1));
set(handles.slider1m1,'Value',s.xc(2));
set(handles.slider10,'Value',s.xc(3));
set(handles.slider11,'Value',s.xc(4));
set(handles.slider2m2,'Value',s.xc(5));
set(handles.slider2m1,'Value',s.xc(6));
set(handles.slider20,'Value',s.xc(7));
set(handles.slider21,'Value',s.xc(8));
set(handles.slider22,'Value',s.xc(9));
set(handles.slider3m3,'Value',s.xc(10));
set(handles.slider3m2,'Value',s.xc(11));
set(handles.slider3m1,'Value',s.xc(12));
set(handles.slider30,'Value',s.xc(13));
set(handles.slider31,'Value',s.xc(14));
set(handles.slider32,'Value',s.xc(15));
set(handles.slider33,'Value',s.xc(16));
set(handles.slider4m4,'Value',s.xc(17));
set(handles.slider4m3,'Value',s.xc(18));
set(handles.slider4m2,'Value',s.xc(19));
set(handles.slider4m1,'Value',s.xc(20));
set(handles.slider40,'Value',s.xc(21));
set(handles.slider41,'Value',s.xc(22));
set(handles.slider42,'Value',s.xc(23));
set(handles.slider43,'Value',s.xc(24));
set(handles.slider44,'Value',s.xc(25));
set(handles.slider5m5,'Value',s.xc(26));
set(handles.slider5m4,'Value',s.xc(27));
set(handles.slider5m3,'Value',s.xc(28));
set(handles.slider5m2,'Value',s.xc(29));
set(handles.slider5m1,'Value',s.xc(30));
set(handles.slider50,'Value',s.xc(31));
set(handles.slider51,'Value',s.xc(32));
set(handles.slider52,'Value',s.xc(33));
set(handles.slider53,'Value',s.xc(34));
set(handles.slider54,'Value',s.xc(35));
set(handles.slider55,'Value',s.xc(36));
set(handles.slider6m6,'Value',s.xc(37));
set(handles.slider6m5,'Value',s.xc(38));
set(handles.slider6m4,'Value',s.xc(39));
set(handles.slider6m3,'Value',s.xc(40));
set(handles.slider6m2,'Value',s.xc(41));
set(handles.slider6m1,'Value',s.xc(42));
set(handles.slider60,'Value',s.xc(43));
set(handles.slider61,'Value',s.xc(44));
set(handles.slider62,'Value',s.xc(45));
set(handles.slider63,'Value',s.xc(46));
set(handles.slider64,'Value',s.xc(47));
set(handles.slider65,'Value',s.xc(48));
set(handles.slider66,'Value',s.xc(49));

guidata(hObject,handles);


% --- Executes on button press in pushbutton_open_figure.
function pushbutton_open_figure_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_open_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
obj = handles.d;

%% generate the shape and scalar field
X = handles.Y_LK(:,1:length(handles.d.xc))* [handles.d.xc(:) handles.d.yc(:) handles.d.zc(:)];
sf = handles.Y_LK(:,1:length(handles.s.xc))*handles.s.xc(:);
%% do the actual plotting
figure;patch('Vertices', X, 'Faces', handles.C,'FaceVertexCData',sf,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);daspect([1 1 1]);axis off;lighting gouraud;


% --- Executes on button press in pushbutton_import_shape.
function pushbutton_import_shape_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_import_shape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiopen('import shape');
handles.d.X_o = X_o;
handles.d = update(handles.d);
plot_dros(handles);
guidata(hObject,handles);