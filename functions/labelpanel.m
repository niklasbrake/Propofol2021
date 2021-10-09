function A = labelpanel(x,y,str)


% str = lower(str); % NPG
str = upper(str); % CellPress

A = annotation('textbox', [x,y,0.05,0.05],'String',str, 'LineStyle', ...
	'none', 'FontWeight','bold', 'FontSize',8,'Margin',0);