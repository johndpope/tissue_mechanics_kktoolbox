function write_ultrascope_xml(fn,S)
% make a basic xml file for the ultrascope configuration
% write the struct s as a basic xml file.
% Uses: struct2xml (included from Matlab file exchange)
% Author: Khaled Khairy (Keller lab)

s.push_config.Attributes.version = '1.0';
s.push_config.info{1}.Attributes.data_root = S.data_root;
s.push_config.info{2}.Attributes.output_root = S.output_root;
s.push_config.info{3}.Attributes.dimensions = S.dimensions;
s.push_config.info{3}.Text = 'in voxel units';
s.push_config.info{4}.Attributes.time_point = S.time_point;
s.push_config.info{5}.Attributes.angle = S.angle;
s.push_config.info{5}.Text = 'in degrees: 0 for c1 and c2, 90 for c1 and c2 when sample rotated';
s.push_config.info{6}.Attributes.channel = S.channel;
s.push_config.info{7}.Attributes.detection_objective = S.detection_objective;
s.push_config.info{8}.Attributes.structured_illumination = S.structured_illumination;
s.push_config.info{9}.Attributes.z_step = S.z_step;
s.push_config.info{9}.Text = 'z step given in microns';
s.push_config.action{1}.Attributes.registration = S.registration;
s.push_config.action{2}.Attributes.fusion = S.fusion;
s.push_config.action{3}.Attributes.deconvolution = S.deconvolution;
s.push_config.action{4}.Attributes.si_reconstruction = S.si_reconstruction;
s.push_config.action{5}.Attributes.segmentation = S.segmentation;
s.push_config.action{6}.Attributes.tracking = S.tracking;
s.push_config.action{7}.Attributes.image_compression = S.image_compression;
s.push_config.action{8}.Attributes.object_space_compression = S.object_space_compression;

xml = struct2xml(s);fid = fopen(fn,'w');fwrite(fid,xml);fclose(fid);
end
%%
function varargout = struct2xml( s, varargin )
%Convert a MATLAB structure into a xml file 
% [ ] = struct2xml( s, file )
% xml = struct2xml( s )
%
% A structure containing:
% s.XMLname.Attributes.attrib1 = 'Some value';
% s.XMLname.Element.Text = 'Some text';
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = '2';
% s.XMLname.DifferentElement{1}.Text = 'Some more text';
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = '2';
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = '1';
% s.XMLname.DifferentElement{2}.Text = 'Even more text';
%
% Will produce:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Written by W. Falkena, ASTI, TUDelft, 27-08-2010
% On-screen output functionality added by P. Orth, 01-12-2010
    
    if (nargin ~= 2)
        if(nargout ~= 1 || nargin ~= 1)
            error(['Supported function calls:' sprintf('\n')...
                   '[ ] = struct2xml( s, file )' sprintf('\n')...
                   'xml = struct2xml( s )']);
        end
    end

    if(nargin == 2)
        file = varargin{1};

        if (isempty(file))
            error('Filename can not be empty');
        end

        if (isempty(strfind(file,'.xml')))
            file = [file '.xml'];
        end
    end
    
    if (~isstruct(s))
        error([inputname(s) ' is not a structure']);
    end
    
    if (length(fieldnames(s)) > 1)
        error(['Error processing the structure:' sprintf('\n') 'There should be a single field in the main structure.']);
    end
    xmlname = fieldnames(s);
    xmlname = xmlname{1};

    %create xml structure
    docNode = com.mathworks.xml.XMLUtils.createDocument(xmlname);

    %process the rootnode
    docRootNode = docNode.getDocumentElement;

    %append childs
    parseStruct(s.(xmlname),docNode,docRootNode,[inputname(1) '.' xmlname '.']);

    if(nargout == 0)
        %save xml file
        xmlwrite(file,docNode);
    else
        varargout{1} = xmlwrite(docNode);
    end  
end

% ----- Subfunction parseStruct -----
function [] = parseStruct(s,docNode,curNode,pName)
    
    fnames = fieldnames(s);
    for i = 1:length(fnames)
        curfield = fnames{i};
        
        if (strcmp(curfield,'Attributes'))
            %Attribute data
            if (isstruct(s.(curfield)))
                attr_names = fieldnames(s.Attributes);
                for a = 1:length(attr_names)
                    cur_attr = attr_names{a};
                    [cur_str,succes] = val2str(s.Attributes.(cur_attr));
                    if (succes)
                        curNode.setAttribute(cur_attr,cur_str);
                    else
                        disp(['Warning. The text in ' pName curfield '.' cur_attr ' could not be processed.']);
                    end
                end
            else
                disp(['Warning. The attributes in ' pName curfield ' could not be processed.']);
                disp(['The correct syntax is: ' pName curfield '.attribute_name = ''Some text''.']);
            end
        elseif (strcmp(curfield,'Text'))
            %Text data
            [txt,succes] = val2str(s.Text);
            if (succes)
                curNode.appendChild(docNode.createTextNode(txt));
            else
                disp(['Warning. The text in ' pName curfield ' could not be processed.']);
            end
        else
            %Sub-element
            if (isstruct(s.(curfield)))
                %single element
                curElement = docNode.createElement(curfield);
                curNode.appendChild(curElement);
                parseStruct(s.(curfield),docNode,curElement,[pName curfield '.'])
            elseif (iscell(s.(curfield)))
                %multiple elements
                for c = 1:length(s.(curfield))
                    curElement = docNode.createElement(curfield);
                    curNode.appendChild(curElement);
                    if (isstruct(s.(curfield){c}))
                        parseStruct(s.(curfield){c},docNode,curElement,[pName curfield '{' num2str(c) '}.'])
                    else
                        disp(['Warning. The cell ' pName curfield '{' num2str(c) '} could not be processed, since it contains no structure.']);
                    end
                end
            else
                %eventhough the fieldname is not text, the field could
                %contain text. Create a new element and use this text
                curElement = docNode.createElement(curfield);
                curNode.appendChild(curElement);
                [txt,succes] = val2str(s.(curfield));
                if (succes)
                    curElement.appendChild(docNode.createTextNode(txt));
                else
                    disp(['Warning. The text in ' pName curfield ' could not be processed.']);
                end
            end
        end
    end
end

% ----- Subfunction val2str -----
function [str,succes] = val2str(val)
    
    succes = true;
    str = [];
    
    if (isempty(val))
        %do nothing
    elseif (ischar(val))
        for i = 1:size(val,1) %multiline string
            str = [str regexprep(val(i,:),'[ ]*', ' ') sprintf('\n')];
        end
        str = str(1:end-1); %skip last enter
    elseif (isnumeric(val))
        tmp = num2str(val);
        for i = 1:size(tmp,1) %multiline string
            str = [str regexprep(tmp(i,:),'[ ]*', ' ') sprintf('\n')];
        end
        str = str(1:end-1); %skip last enter
    else
        succes = false;        
    end
end