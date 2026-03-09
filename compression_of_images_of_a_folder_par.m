% Set the directory path
image_dir = '../URBAN100';

% Get all image files (adjust extensions as needed)
extensions = {'*.png', '*.jpg', '*.tiff', '*.bmp'};
files = [];
for i = 1:length(extensions)
    files = [files; dir(fullfile(image_dir, extensions{i}))];
end

%% Loop over all images
images = files(1:end);
nImages = numel(images);
nReps   = 1;

% Pre-allocate the full results cell array
RESULTS_GLOBAL = cell(nImages, nReps);

parfor i = 1:nImages
    filename = images(i).name;
    filepath = fullfile(image_dir, filename);
    disp(filepath)

    % Run compression nReps times for this image
    results_i = cell(1, nReps);
    for k = 1:nReps
        results_i{k} = imagecompression_color_OW_par(filepath);
    end
    RESULTS_GLOBAL(i, :) = results_i;
end

% save('RANDOM_OW_par2.mat', 'RESULTS_GLOBAL');