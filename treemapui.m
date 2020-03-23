function rectangles = treemapui(data,w,h)
    %TREEMAP   Partition a rectangular area.
    %   rectangles = treemap(data,w,h) partitions the rectangle [0,0,w,h]
    %   into an array of rectangles, one for each element of the vector "data"
    %   The areas of the rectangles are proportional to the values in data and
    %   the aspect ratio of the rectangles are close to one.
    %
    %   The algorithm used is as follows. The data are sorted from largest to
    %   smallest. We then lay out a row or column along the shortest side (w or h)
    %   of the remaining rectangle. Blocks are added to this new row or column
    %   until adding a new block would worsen the average aspect ratio of the
    %   blocks in the row or column. Once this row or column is laid out, we
    %   recurse on the remaining blocks and the remaining rectangle.
    %
    %   Examples:
    %    treemap(rand(20,1))
    %
    %    r = treemap(1:15,1.6,1);
    %    plotrectangles(r)
    %    outline(r)
    %
    %   Copyright 2007-2013 The MathWorks, Inc.
    %

    % default recursion limit of 500 is for wimps.
    set(0,'RecursionLimit',5000)

    if (nargin == 0)
        % make up some random data
        data = abs(randn(1,80)) .^ 2;
    end

    % if you don't specify a rectangle, we use the unit square
    if(nargin < 2)
        w = 1;
        h = 1;
    end

    % this is the ratio of rectangle area to data values
    areaRatio = (w * h) / sum(data);

    % we lay out the date from largest to smallest
    [v,sortIndex] = sort(data,'descend');

    % this chooses the number of blocks in each row or column
    p = partitionArea(v,w,h);

    % this leaves us with the task of assigning rectangles to each element.

    % storage for the result
    rectangles = zeros(4,length(data));

    % these hold the origin for each row
    oX = 0;
    oY = 0;

    % the index of the first value in v for the current row
    first = 1;
    % where to put the resulting rectangles in the rectangles array
    resultIndex = 1;
    % for each row to layout...
    for i = 1:length(p)
        % which values are we placing
        last = first + p(i) - 1;
        chunk = v(first:last);
        % for the next iteration
        first = last + 1;

        % how much area should each block have?
        blockArea = chunk * areaRatio;

        % how much area must the entire column have?
        columnArea = sum(blockArea);

        % the origin for the blocks starts as the origin for the column
        blockX = oX;
        blockY = oY;

        % are we laying out a row or a column?
        if((w < h)) % a row
            % so how thick must the row be?
            columnHeight = columnArea / w;
            % lets place each value
            for j = 1:p(i)
                % so how wide should it be?
                blockWidth = blockArea(j) / columnHeight;
                % remember the rectangle
                rectangles(:,sortIndex(resultIndex)) = [blockX,blockY,blockWidth,columnHeight];
                resultIndex = resultIndex + 1;
                % move the origin for the nextblock
                blockX = blockX + blockWidth;
            end
            % move the origin for the next row
            oY = oY + columnHeight;
            h = h - columnHeight;
        else % a column
            columnWidth = columnArea / h;
            % lets place each value
            for j = 1:p(i)
                % so how high should it be?
                blockHeight = blockArea(j) / columnWidth;
                % if one corner is at oX,oY, where is the opposite corner?
                rectangles(:,sortIndex(resultIndex)) = [blockX,blockY,columnWidth,blockHeight];
                resultIndex = resultIndex + 1;
                % move the origin for the nextblock
                blockY = blockY  +  blockHeight;
            end
            % move the origin for the next row
            oX = oX + columnWidth;
            w = w - columnWidth;
        end
    end

    % show the results if they are not returned
    if(nargout == 0)
        plotRectanglesui(rectangles)
    end
end

function partition = partitionArea(v,w,h)
    % return a vector that tells us how many blocks go in each column or
    % row. Sum(partition) == length(v);

    % the rest of the code only wands to deal with long side and short
    % side. It is not interested in orientation (w or h)
    if(w>h)
        longSide = w;
        shortSide = h;
    else
        longSide = h;
        shortSide = w;
    end

    % this is the ratio of value units to associated area.
    areaRatio = (w * h) / sum(v);

    % we want to minimize cost
    bestAverageAspectRatio = inf;

    % we will return an array of integers that tell how many blocks to
    % place in each row (or column)
    partition = [];

    % How many blocks should go in the next column of rectangles?
    % i is our current guess. We'll keep incrementing it until aspect ratio
    % (cost) gets worse.
    for i = 1:length(v)
        columnTotal = sum(v(1:i));
        columnArea = columnTotal * areaRatio;
        columnWidth = columnArea / shortSide;

        % this is the cost associated with this value of i.
        sumOfAspectRatio = 0;
        % for each block in the column
        for j = 1:i
            % what is the aspect ratio
            blockArea = v(j) * areaRatio;
            blockHeight = blockArea / columnWidth;
            aspectRatio = blockHeight / columnWidth;
            if(aspectRatio < 1)
                aspectRatio = 1 / aspectRatio;
            end
            sumOfAspectRatio = sumOfAspectRatio + aspectRatio;
        end

        averageAspectRatio = sumOfAspectRatio/i;

        % if the cost of this value of i worse than for (i-1) we are done
        % and we will use i-i blocks in this column.
        if(averageAspectRatio >= bestAverageAspectRatio)
            if(partition < i) % if we are not done yet, we recurse
                p = partitionArea(v(i:end),shortSide,longSide - columnWidth);
                partition = [partition,p];
            end
            return
        end

        bestAverageAspectRatio = averageAspectRatio;
        partition = i;
    end

end
