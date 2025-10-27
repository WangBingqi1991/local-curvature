analyze_vessel6("your.tif")
function analyze_vessel6(image_path)
    if isstring(image_path)
        image_path = char(image_path);
    end
    
    if ~exist(image_path, 'file')
        error('File not found: %s', image_path);
    end
    img = imread(image_path);
    
    figure;
    imshow(img);
    
    img_hsv = rgb2hsv(img);
    
    hue_red1 = img_hsv(:,:,1) > 0.95;
    hue_red2 = img_hsv(:,:,1) < 0.05;
    saturation = img_hsv(:,:,2) > 0.5;
    value = img_hsv(:,:,3) > 0.3;
    
    red_mask = (hue_red1 | hue_red2) & saturation & value;
    
    red_channel = img(:,:,1);
    green_channel = img(:,:,2);
    blue_channel = img(:,:,3);
    
    red_diff_mask = (red_channel - green_channel > 50) & ...
                   (red_channel - blue_channel > 50) & ...
                   (red_channel > 100);
    
    combined_mask = red_mask | red_diff_mask;
    
    figure;
    imshow(combined_mask);
    
    clean_mask = bwareaopen(combined_mask, 100);
    clean_mask = imclose(clean_mask, strel('disk', 5));
    
    edge_mask = bwperim(clean_mask, 8);
    
    if nnz(edge_mask) < 50
        gray_img = rgb2gray(img);
        edge_mask = edge(gray_img, 'Canny', [0.1, 0.3]);
    end
    
    [edge_y, edge_x] = find(edge_mask);
    
    if isempty(edge_x) || numel(edge_x) < 10
        error('Insufficient edge points detected');
    end
    
    figure;
    imshow(img); hold on;
    plot(edge_x, edge_y, 'r.', 'MarkerSize', 5);
    
    [xc, yc, R] = circfit(edge_x, edge_y);
    
    theta_deg = (0:359)';
    theta_rad = deg2rad(theta_deg);
    
    avg_x = xc + R * cos(theta_rad);
    avg_y = yc + R * sin(theta_rad);
    
    dx = edge_x - xc;
    dy = edge_y - yc;
    dist_real = sqrt(dx.^2 + dy.^2);
    theta_real = atan2(dy, dx);
    
    theta_real_deg = mod(rad2deg(theta_real), 360);
    
    dist_diff = zeros(360, 1);
    bin_counts = zeros(360, 1);
    
    for i = 1:length(edge_x)
        bin_idx = floor(theta_real_deg(i)) + 1;
        if bin_idx > 360
            bin_idx = 1;
        end
        dist_diff(bin_idx) = dist_diff(bin_idx) + (dist_real(i) - R);
        bin_counts(bin_idx) = bin_counts(bin_idx) + 1;
    end
    
    for k = 1:360
        if bin_counts(k) > 0
            dist_diff(k) = dist_diff(k) / bin_counts(k);
        else
            prev_idx = mod(k-2, 360) + 1;
            next_idx = mod(k, 360) + 1;
            dist_diff(k) = (dist_diff(prev_idx) + dist_diff(next_idx)) / 2;
        end
    end
    
    result_table = table(theta_deg, dist_diff, ...
        'VariableNames', {'Angle_Deg', 'Distance_Difference'});
    
    [folder, name, ~] = fileparts(image_path);
    if isempty(folder)
        folder = pwd;
    end
    excel_path = fullfile(folder, [name '_analysis.xlsx']);
    
    writetable(result_table, excel_path);
    
    figure;
    imshow(img); hold on;
    plot(avg_x, avg_y, 'g-', 'LineWidth', 2);
    plot(edge_x, edge_y, 'r.', 'MarkerSize', 8);
    plot(xc, yc, 'b+', 'MarkerSize', 20, 'LineWidth', 3);
    
    figure;
    plot(theta_deg, dist_diff, 'b-', 'LineWidth', 1.5);
    xlabel('Angle (deg)');
    ylabel('Distance Difference (px)');
    grid on;
    
    disp(['Results saved to: ' excel_path]);
end

function [xc, yc, R] = circfit(x, y)
    x = x(:); y = y(:);
    
    mx = mean(x);
    my = mean(y);
    x_centered = x - mx;
    y_centered = y - my;
    
    Z = [x_centered, y_centered, ones(size(x))];
    H = [sum(x_centered.^2), sum(x_centered.*y_centered), sum(x_centered);
         sum(x_centered.*y_centered), sum(y_centered.^2), sum(y_centered);
         sum(x_centered), sum(y_centered), numel(x)];
    B = -[sum(x_centered.*(x_centered.^2 + y_centered.^2));
          sum(y_centered.*(x_centered.^2 + y_centered.^2));
          sum(x_centered.^2 + y_centered.^2)];
    
    a = H \ B;
    
    xc = -a(1)/2 + mx;
    yc = -a(2)/2 + my;
    R = sqrt(a(1)^2/4 + a(2)^2/4 - a(3));
end