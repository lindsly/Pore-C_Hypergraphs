hs_cytobands = cytobandread('hs_cytoBand.txt');
% hfchrom = chromosomeplot_sl(hs_cytobands);
hfchrom = chromosomeplot(hs_cytobands);

set(gcf, 'color', 'none'); set(gca, 'color', 'none');
fig = get(gca);

%save as rendered output with transparent background not work
exportgraphics(fig, 'test.png', 'ContentType', 'image', 'BackGroundColor', 'none');
% prints warning: Warning: Background transparency is not supported; using white instead.

% now our solution, export with two background colors
exportgraphics(fig, 'test1.png', 'ContentType', 'image', 'BackgroundColor', 'k'); % black background
exportgraphics(fig, 'test2.png', 'ContentType', 'image', 'BackgroundColor', 'w'); % white background

% load exported images back in and scale to [0,1]
u = imread('test1.png');
u = double(u) / 255;
v = imread('test2.png');
v = double(v) / 255;

% recover transparency as a
a = 1 - v + u;
a = mean(a, 3);
a = max(0, min(1, a));
m = a > eps;

% recover rgb
c = zeros(size(u));
for i = 1 : 3
    ui = u(:, :, i);
    ci = c(:, :, i);
    ci(m) = ui(m) ./ a(m);
    c(:, :, i) = ci;
end
c = max(0, min(1, c));

% store again
imwrite(uint8(c*255), 'test.transparent.png', 'Alpha', a);

% temporary files test1.png and test2.png can now be deleted