function handle = findFigure(figName)
    figHandles = get(groot, 'Children');
    
    figExist = false;
    for i = 1:length(figHandles)
        if isequal(figHandles(i).Name, figName)
            figExist = true;
            handle = figure(figHandles(i));
            break;
        end
    end
    
    if ~figExist
        handle = figure('Name', figName);
    end
end