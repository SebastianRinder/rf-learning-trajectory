function skipVisValue = selectFigure(figName)
    figHandles = get(groot, 'Children');
    
    figExist = false;
    for i = 1:length(figHandles)
        if isequal(figHandles(i).Name, figName)
            figExist = true;
            myFigure = figure(figHandles(i));
            c = findobj(myFigure, 'type', 'uicontrol', 'style', 'checkbox');    
            skipVisValue = get(c, 'Value');
            break;
        end
    end
    
    if ~figExist
        figure('Name', figName);
        skipVisValue = 0;
    end
end