% Step 1: Generate synthetic uncompressed grayscale video
v = VideoWriter('test_raw.avi', 'Uncompressed AVI');
v.FrameRate = 25;
open(v);

nFrames = 30;
H = 240;
W = 320;

for i = 1:nFrames
    % Create a moving gradient as a simple synthetic video
    [X, Y] = meshgrid(1:W, 1:H);
    frame = uint8(mod(X + Y + i*5, 256));  % shifting pattern
    frame_rgb = cat(3, frame, frame, frame); % convert to RGB (required by VideoWriter)
    writeVideo(v, frame_rgb);
end
close(v);
disp('Video written successfully')

% Step 2: Read it back as a 3D tensor
vid = VideoReader('test_raw.avi');

nFrames = ceil(vid.Duration * vid.FrameRate) + 1;
tensor = zeros(H, W, nFrames, 'uint8');

i = 1;
while hasFrame(vid)
    frame = readFrame(vid);
    if size(frame, 3) == 3
        frame = rgb2gray(frame);
    end
    tensor(:,:,i) = frame;
    i = i + 1;
end
tensor = tensor(:,:,1:i-1);  % trim unused frames

fprintf('Tensor size: %d x %d x %d\n', size(tensor,1), size(tensor,2), size(tensor,3));

% Step 3: Visualize a few frames
figure;
for k = 1:6
    subplot(2,3,k);
    imshow(tensor(:,:,k));
    title(['Frame ' num2str(k)]);
end
% ```
% 
% ## What This Does
% 
% - **Generates** a 30-frame synthetic video with a shifting gradient pattern — simple but shows motion
% - **Saves** it as a truly uncompressed AVI using MATLAB's `Uncompressed AVI` profile
% - **Reads** it back into a 3D tensor of shape **(H x W x T)**
% - **Visualizes** the first 6 frames in a figure
% 
% ## Expected Output
% ```
% Video written successfully
% Tensor size: 240 x 320 x 30