function [BW] = make_circle (N, M, x_center, y_center, the_radius)
    % make_circle(N,M, the_radius):
    %   Creates an ellipse of size NxM.
    [x, y]       = meshgrid(1:M, 1:N);  % Reverse order.
    BW = (x - x_center).^2 + (y - y_center).^2 <= the_radius^2;
end
