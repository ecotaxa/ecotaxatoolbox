function plotRectangles(rectangles,labels,colors)
    %PLOTRECTANGLES  Create patches to display rectangles
    %   plotRectangles(rectangles,labels,colors)
    %   rectangles is a 4-by-N matrix where each column is a rectangle in
    %   the form [x0,y0,w,h]. x0 and y0 are the coordinates of the lower
    %   left-most point in the rectangle.
    %
    %   Example:
    %    cla
    %    plotRectangles([0 0 1 1; 1 0 1 1; 0 1 2 1]')
    %
    %   Copyright 2007-2013 The MathWorks, Inc.
    %
    %   See also OUTLINE.
    
    if(nargin < 2)
        labels = [];
    end
    
    if(nargin < 3)
        colors = rand(size(rectangles,2),3).^0.5;
    end
    
    % make a patch for each rectangle
    
    for i = 1:size(rectangles,2)
        r = rectangles(:,i);
        xPoints = [r(1), r(1), r(1) + r(3),r(1) + r(3)];
        yPoints = [r(2), r(2)+ r(4), r(2)+ r(4),r(2)];
        patch(xPoints,yPoints,colors(i,:),'EdgeColor','none');
        if(~isempty(labels))
            text(r(1) + r(3)/2,r(2) + r(4)/2, 1, labels{i}, ...
                'VerticalAlignment','middle','HorizontalAlignment','center')
        end
    end
    
    axis equal
    axis tight
    axis off
    
end