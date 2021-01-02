function idx_choice = ROIClassChooseDlg(ROIstring,varargin)

if isempty(varargin)
    XYPos = [0.5, 0.5];
else
    XYPos = varargin{1};
end
Pos = [XYPos, 0.25,0.1];
Pos(1) = Pos(1) - (Pos(3)/2);
Pos(2) = Pos(2) - (Pos(4)/2);
d = dialog('Units','normalized','Position',Pos,'Name', 'Select one');

txt = uicontrol('Parent',d,'Style','text','Units','normalized','Position',[0.25, 0.35, 0.5, 0.5],...
    'String','Choose an Class to apply to the current ROI');

popup = uicontrol('Parent',d, 'Style','popup','Units','normalized','Position', [0.25 0.25 0.5 0.5], ...
    'String',ROIstring,'Callback',@ROIchoicePop_callback);

btn = uicontrol('Parent',d,'Units','normalized','Position',[0.35,0.1,0.25,0.25], 'String', 'Select','Callback','delete(gcf)');

idx_choice = 1;
uiwait(d);

    function ROIchoicePop_callback(popup,callbackdata)
        
        idx_choice = get(popup,'Value');
        
    end
end