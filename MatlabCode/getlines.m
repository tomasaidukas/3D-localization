function [x, y] = getlines(centre, radius, ang)

% Define x and y end co-ordinates
endx = cos(ang) * radius + centre;
endy = cos(pi / 2 - ang) * radius + centre;

% Generate a line along an angle
% y = mx + b
% b = centre
slope = (endy - centre) / (endx - centre);
intercept = centre * (1 - slope);

x = linspace(centre, endx, radius);
y = slope .* x + intercept;

if ang == pi/2
    x = centre
    y = linspace(centre, endy, radius);
end

end