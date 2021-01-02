function idx_choice = ROIchooseDlg(ROIstring)

XYPos = [0.5, 0.5];
Pos = [XYPos, 0.1,0.075];
Pos(1) = Pos(1) - (Pos(3)/2);
Pos(2) = Pos(2) - (Pos(4)/2);
d = dialog('Units','normalized','Position',Pos,'Name', 'Select one');

txt = uicontrol('Parent',d,'Style','text','Units','normalized','Position',[0.05 0.60 0.9 0.33],...
    'String','Choose an ROI to Analyze');

popup = uicontrol('Parent',d, 'Style','popup','Units','normalized','Position', [0.3 0.58 0.4 0.20], ...
    'String',ROIstring,'Callback',@ROIchoicePop_callback);

btn = uicontrol('Parent',d,'Units','normalized','Position',[0.35 0.17 0.27 0.2], 'String', 'Select','Callback','delete(gcf)');

idx_choice = 1;
uiwait(d);

    function ROIchoicePop_callback(popup,callbackdata)
        
        idx_choice = get(popup,'Value');
        
    end
end