function S = read_ultrascope_xml(fn)
% Reads ultrascope data acquisition information and instructions on how to
% proceed from an xml file. Returns a matlab struct.
% Author: Khaled Khairy (Keller lab)
%% parse the XML file and reconstruct the input struct
S = struct( 'data_header',[],...
            'output_root',[],...
            'specimen_name',[], ...
            'time_point',[], ...
            'angle', [],...
            'camera_index', [],...
            'wavelength', [], ...
            'illumination_filter', [],...
            'detection_filter', [], ...
            'dimensions',[],...
            'z_step', [], ...
            'planes', [], ...
            'detection_objective', [],...
            'SI', [], ...
            'experiment_notes', [], ...
            'custom_pre_script', [], ...
            'custom_pre_commands', [], ...
            'custom_post_script', [], ...
            'custom_post_commands', [], ...
            'registration', [], 'fusion', [], 'deconvolution', [], ...
            'si_reconstruction', [], 'segmentation', [], 'tracking', [], ...
            'image_compression', [], 'object_space_compression', [],'v3d',[]);

SS = parseXML(fn);  % read in the xml file
cc = SS.Children;
%% parse the information and construct the final matlab struct
info_fields = fieldnames(S);
for ix = 1:numel(info_fields)   % loop over the fields we want to find
    for jx = 1:numel(cc)    % loop over the elements
        if ~isempty(cc(jx).Attributes),
            if strcmp(cc(jx).Attributes.Name, info_fields{ix})
                S = setfield(S, info_fields{ix}, cc(jx).Attributes.Value);
            end
        end
    end
end
% % %%% convert datatypes
% % S.dimensions = str2num(S.dimensions);
% % S.angle = str2num(S.angle);
% % S.z_step = str2num(S.z_step);


%%
function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
    tree = xmlread(filename);
catch
    error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems
% with very deeply nested trees.
try
    theStruct = parseChildNodes(tree);
catch
    error('Unable to parse XML file %s.',filename);
end


% ----- Subfunction PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
    childNodes = theNode.getChildNodes;
    numChildNodes = childNodes.getLength;
    allocCell = cell(1, numChildNodes);
    
    children = struct(             ...
        'Name', allocCell, 'Attributes', allocCell,    ...
        'Data', allocCell, 'Children', allocCell);
    
    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end

% ----- Subfunction MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
    'Name', char(theNode.getNodeName),       ...
    'Attributes', parseAttributes(theNode),  ...
    'Data', '',                              ...
    'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
    nodeStruct.Data = char(theNode.getData);
else
    nodeStruct.Data = '';
end

% ----- Subfunction PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
    theAttributes = theNode.getAttributes;
    numAttributes = theAttributes.getLength;
    allocCell = cell(1, numAttributes);
    attributes = struct('Name', allocCell, 'Value', ...
        allocCell);
    
    for count = 1:numAttributes
        attrib = theAttributes.item(count-1);
        attributes(count).Name = char(attrib.getName);
        attributes(count).Value = char(attrib.getValue);
    end
end
