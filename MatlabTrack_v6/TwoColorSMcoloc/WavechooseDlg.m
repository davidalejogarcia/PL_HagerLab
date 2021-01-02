function idx_choice = WavechooseDlg(wstring)

d = dialog('Position',[500 500 250 160],'Name', 'Select one');

txt = uicontrol('Parent',d,'Style','text','Position',[20 90 210 40],...
    'String','Select the wavelength of each channel');
txt1 = uicontrol('Parent',d,'Style','text','Position',[-25 50 210 40],...
    'String','Channel 1:');
txt2 = uicontrol('Parent',d,'Style','text','Position',[-25 10 210 40],...
    'String','Channel 2:');

popup1 = uicontrol('Parent',d, 'Style','popup','Position', [115 70 70 25], ...
    'String',wstring,'Callback',@WchoicePop1_callback);

popup2 = uicontrol('Parent',d, 'Style','popup','Position', [115 30 70 25], ...
    'String',wstring,'Callback',@WchoicePop2_callback);


btn = uicontrol('Parent',d,'Position',[89 5 70 25], 'String', 'Select','Callback','delete(gcf)');

idx_choice = [1 1];
uiwait(d);

    function WchoicePop1_callback(popup,callbackdata)
        
        idx_choice(1) = get(popup,'Value');
        
    end

function WchoicePop2_callback(popup,callbackdata)
        
        idx_choice(2) = get(popup,'Value');
        
    end
end