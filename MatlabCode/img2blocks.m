function image3D = img2blocks(img, rowsize, colsize)

% Takes an image and splits it into blocks

    [rows columns numberOfColorBands] = size(img);
    % Let's assume we know the block size and that all blocks will be the same size.
    blockSizeR = rowsize; % Rows in block.
    blockSizeC = colsize; % Columns in block.

    % Figure out the size of each block. 
    wholeBlockRows = floor(rows / blockSizeR);
    wholeBlockCols = floor(columns / blockSizeC);

    % Preallocate a 3D image
    image3d = zeros(wholeBlockRows, wholeBlockCols, 36);

    % Now scan though, getting each block and putting it as a slice of a 3D array.
    sliceNumber = 1;
    for row = 1 : blockSizeR : rows
        for col = 1 : blockSizeC : columns
            % Let's be a little explicit here in our variables
            % to make it easier to see what's going on.
            row1 = row;
            row2 = row1 + blockSizeR - 1;
            col1 = col;
            col2 = col1 + blockSizeC - 1;
            % Extract out the block into a single subimage.
            oneBlock = img(row1:row2, col1:col2);

            % Assign this slice to the image we just extracted.
            image3D(:, :, sliceNumber) = oneBlock;
            sliceNumber = sliceNumber + 1;
        end
    end
end