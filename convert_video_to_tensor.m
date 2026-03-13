v = VideoReader('movies/8052353-hd_1920_1080_30fps.mp4');

H = v.Height;
W = v.Width;
nFrames = floor(v.Duration * v.FrameRate);

tensor = zeros(H, W, nFrames, 'uint8');

i = 1;
while hasFrame(v)
    frame = readFrame(v);
    if size(frame, 3) == 3
        frame = rgb2gray(frame);
    end
    tensor(:,:,i) = frame;
    i = i + 1;
end

disp(size(tensor))  % (H, W, T)

I = tensor;

save('movies/movie03.mat','I')