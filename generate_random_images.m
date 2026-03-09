function generate_random_images(max_power)
    % Create output directory if it doesn't exist
    output_dir = '../random_images';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Loop through powers of 2 from 2^1 to 2^max_power
    for p = 4:max_power
        size_n = 2^p;
        
        % Generate random RGB image (height x width x 3 channels), uint8 [0,255]
        img = uint8(randi([0, 255], size_n, size_n, 3));
        
        % Build filename, e.g. "image_2^3_8x8.png"
        filename = fullfile(output_dir, sprintf('image_2pow%d_%dx%d.png', p, size_n, size_n));
        
        % Write image to file
        imwrite(img, filename);
        
        fprintf('Saved: %s\n', filename);
    end
    
    fprintf('Done! %d images saved to %s\n', max_power, output_dir);
end