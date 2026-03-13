function [] = write_tensor_into_video(tensor,filename)
% Assume tensor is (H x W x T) of uint8
v = VideoWriter(filename, 'Uncompressed AVI');
v.FrameRate = 30;
open(v);

for i = 1:size(tensor, 3)
    frame = uint8(tensor(:,:,i));                  % extract single grayscale frame (H x W)
    frame_rgb = cat(3, frame, frame, frame); % replicate to RGB (required by VideoWriter)
    writeVideo(v, frame_rgb);
end

close(v);
disp('Video saved successfully')