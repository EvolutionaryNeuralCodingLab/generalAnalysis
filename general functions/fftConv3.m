function out = fftConv3(videoTrials, kernel)
%FFTCONV3 Convolve video along 3rd dimension with 1D kernel using FFT
%
%   out = fftConv3(videoTrials, kernel)
%
%   videoTrials : [Y x X x T] array
%   kernel      : [1 x 1 x K] or [K] kernel along time
%   out         : [Y x X x T] convolution result, 'same' size

    % Ensure kernel is column vector along 3rd dimension
    if isvector(kernel)
        kernel = reshape(kernel,1,1,[]);
    elseif size(kernel,1) ~= 1 || size(kernel,2) ~= 1
        error('Kernel must be 1x1xK or vector.');
    end

    vidSize = size(videoTrials,3);
    kerSize = size(kernel,3);
    fftLen  = vidSize + kerSize - 1;

    % FFT along 3rd dimension
    Fvid = fft(videoTrials, fftLen, 3);
    Fker = fft(kernel, fftLen, 3);

    % Multiply in frequency domain
    Fprod = Fvid .* Fker;

    % Back to time domain, take real part
    tmp = ifft(Fprod, [], 3, 'symmetric');

    % Crop to 'same'
    out = tmp(:,:,ceil(kerSize/2):ceil(kerSize/2)+vidSize-1);
end
