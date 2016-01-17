function SaveAllFigures(opt,filetype)

if nargin == 0
opt='Unknown';
end
if nargin < 2
filetype = 'fig';
end

ChildList = sort(get(0,'Children'));
for cnum = 1:length(ChildList)
if strncmp(get(ChildList(cnum),'Type'),'figure',6)
saveas(ChildList(cnum), ['FigureSave', opt, '_', num2str(ChildList(cnum)), '.' filetype]);
end
end