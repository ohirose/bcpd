function J = medianfilter(A, mn, padopt)
% MY_MEDFILT2 Mimics standard medfilt2 without the Image Processing Toolbox.
%   J = MY_MEDFILT2(A) uses a 3x3 window and zero padding.
%   J = MY_MEDFILT2(A, [M N]) uses an MxN window.
%   J = MY_MEDFILT2(A, [M N], PADOPT) uses specified padding.
%       PADOPT can be 'zeros' (default) or 'replicate' (repeats border pixels).

    % --- 1. Argument Parsing & Defaults ---
    if nargin < 2
        mn = [3 3]; % Default window size
    end
    if nargin < 3
        padopt = 'zeros'; % Default padding (matches medfilt2 standard)
    end
    
    % Validate input is 2D
    if ~ismatrix(A)
        error('Input must be a 2D matrix.');
    end

    % Store original class to restore later
    inputClass = class(A);
    
    % Work in double precision for calculation
    A = double(A);
    [rows, cols] = size(A);
    
    % Dimensions of the window
    m = mn(1);
    n = mn(2);
    
    % Calculate padding thickness
    pad_r = floor(m / 2);
    pad_c = floor(n / 2);
    
    % --- 2. Manual Padding ---
    % We cannot use padarray (toolbox function), so we do it manually.
    if strcmp(padopt, 'replicate')
        % Replicate borders (slower but usually looks better)
        padded = zeros(rows + 2*pad_r, cols + 2*pad_c);
        % Center
        padded(pad_r+1 : end-pad_r, pad_c+1 : end-pad_c) = A;
        % Top/Bottom
        padded(1:pad_r, pad_c+1:end-pad_c) = repmat(A(1,:), pad_r, 1);
        padded(end-pad_r+1:end, pad_c+1:end-pad_c) = repmat(A(end,:), pad_r, 1);
        % Left/Right (including corners)
        padded(:, 1:pad_c) = repmat(padded(:, pad_c+1), 1, pad_c);
        padded(:, end-pad_c+1:end) = repmat(padded(:, end-pad_c), 1, pad_c);
        
    else 
        % 'zeros' (The default behavior of medfilt2)
        padded = zeros(rows + 2*pad_r, cols + 2*pad_c);
        padded(pad_r+1 : end-pad_r, pad_c+1 : end-pad_c) = A;
    end
    
    % --- 3. Vectorized "Shift and Stack" ---
    % Pre-allocate stack
    % Size: height x width x (window_area)
    stack = zeros(rows, cols, m*n);
    idx = 1;
    
    % Loop through every neighbor in the MxN window
    for r = -pad_r : (m - pad_r - 1) % Handles even/odd windows correctly
        for c = -pad_c : (n - pad_c - 1)
            
            % Calculate window bounds in the padded image
            r_start = 1 + pad_r + r;
            r_end   = rows + pad_r + r;
            c_start = 1 + pad_c + c;
            c_end   = cols + pad_c + c;
            
            % Grab the shifted image layer
            stack(:, :, idx) = padded(r_start:r_end, c_start:c_end);
            idx = idx + 1;
        end
    end
    
    % --- 4. Compute Median & Restore Class ---
    J = median(stack, 3);
    
    % Cast back to original type (e.g., uint8)
    J = cast(J, inputClass);
end
