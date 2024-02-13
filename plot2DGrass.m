function plot2DGrass(Ys, YsNoise, time)

if nargin < 1 || nargin > 3
    error('Please input one or two variables');
end

if nargin == 1
    N = length(Ys);
    if N <= 0
        error('No element to show');
    end
    switch N
        case 1
            line([0 Ys{1}(1)], [0 Ys{1}(2)], 'Color', [1 0 0]);
        case 2
            line([0 Ys{1}(1)], [0 Ys{1}(2)], 'Color', [1 0 0]);
            line([0 Ys{2}(1)], [0 Ys{2}(2)], 'Color', [0 0 1]);
        otherwise
            cmap = jet(N);
            for iI = 1:N
                line([0 Ys{iI}(1)], [0 Ys{iI}(2)], 'Color', cmap(iI, :));
                %subplot(1, 2, 1), line([0 Ys{iI}(1)], [0 Ys{iI}(2)], 'Color', cmap(iI, :)); hold on;
                %subplot(1, 2, 2), scatter(Ys{iI}(1), Ys{iI}(2), 10, cmap(iI, :), 'filled'); hold on;
            end
            %subplot(1, 2, 1), hold off, axis equal;
            %subplot(1, 2, 2), hold off, axis equal;
            axis equal
    end
end

if nargin == 2
    matYs = zeros(length(Ys), 2);
    matYsNoise = zeros(length(YsNoise), 2);
    for iI = 1:length(Ys)
        matYs(iI, 1) = Ys{iI}(1);
        matYs(iI, 2) = Ys{iI}(2);
    end
    for iI = 1:length(YsNoise)
        matYsNoise(iI, 1) = YsNoise{iI}(1);
        matYsNoise(iI, 2) = YsNoise{iI}(2);
    end
    plot(matYs(:, 1), matYs(:, 2), '-k*', 'LineWidth', 3 );
    hold on
    plot(matYsNoise(:, 1), matYsNoise(:, 2), 'rs', 'LineWidth', 3);
    hold off
    axis equal
end

if nargin == 3
    angleYs = zeros(length(Ys), 1);
    angleYsNoise = zeros(length(YsNoise), 1);
    for iI = 1:length(Ys)
        angleYs(iI) = atan(Ys{iI}(2)/Ys{iI}(1));
    end
    for iI = 1:length(YsNoise)
        angleYsNoise(iI) = atan(YsNoise{iI}(2)/YsNoise{iI}(1));
    end
    if ~isempty(angleYs)
    plot(time, angleYs, '-r', 'LineWidth', 3);
    end
    hold on
    if ~isempty(angleYsNoise)
    plot(time, angleYsNoise, 'ks', 'LineWidth', 3);
    end
    hold off
end