%% Comparative Analysis of Volumetric Methods for Deep Excavation in Abuja
% Author: Shekwomagode Aaron Tonya Solomon
% B.Tech. Surveying and Geoinformatics — Federal University of Technology, Minna
%
% Companion MATLAB implementation mirroring `src/volumetric_methods.py`
% Evaluates 5 standard Surveying & Geoinformatics volumetric formulations
% on the Abuja Central Area (UTM Zone 32N) deep excavation dataset.

clear; clc; close all;

%% 1. Load Abuja Survey Grid Dataset
T = readtable('../data/abuja_deep_excavation_survey.csv');

xs = unique(T.Local_X_m)';
ys = unique(T.Local_Y_m)';
nx = numel(xs);
ny = numel(ys);
dx = xs(2) - xs(1);
dy = ys(2) - ys(1);

% Reshape 1D table columns into structured 2D matrices (ny x nx)
XX    = reshape(T.Easting_m,  [nx, ny])';
YY    = reshape(T.Northing_m, [nx, ny])';
Z_pre = reshape(T.Pre_Elev_m, [nx, ny])';
Z_post= reshape(T.Post_Elev_m,[nx, ny])';
D     = reshape(T.Cut_Depth_m,[nx, ny])';

%% 2. Method 1 — Average End-Area (Trapezoidal Cross-Section) Rule
section_areas_trap = trapz(xs, D, 2);
V_end_area = trapz(ys, section_areas_trap);

%% 3. Method 2 — Prismoidal (Composite Simpson's 1/3) Rule
wx = ones(1, nx); wx(2:2:end-1) = 4; wx(3:2:end-2) = 2;
wy = ones(ny, 1); wy(2:2:end-1) = 4; wy(3:2:end-2) = 2;
section_areas_simp = (dx / 3) * (D * wx');
V_prismoidal = (dy / 3) * (wy' * section_areas_simp);

%% 4. Method 3 — Spot-Height Rectangular Grid (1-2-3-4 Borrow-Pit Weights)
W_grid = zeros(ny, nx);
W_grid(1:end-1, 1:end-1) = W_grid(1:end-1, 1:end-1) + 1;
W_grid(1:end-1, 2:end)   = W_grid(1:end-1, 2:end)   + 1;
W_grid(2:end,   1:end-1) = W_grid(2:end,   1:end-1) + 1;
W_grid(2:end,   2:end)   = W_grid(2:end,   2:end)   + 1;
V_grid = ((dx * dy) / 4) * sum(W_grid .* D, 'all');

%% 5. Method 4 — Horizontal Contour Slicing (1.0 m Contour Interval)
ci = 1.0;
levels = 0 : ci : ceil(max(D, [], 'all'));
D_cell = 0.25 * (D(1:end-1,1:end-1) + D(1:end-1,2:end) + D(2:end,1:end-1) + D(2:end,2:end));
A_slice = arrayfun(@(z) sum(D_cell >= z, 'all') * dx * dy, levels);
V_contour = sum((ci / 3) .* (A_slice(1:end-1) + A_slice(2:end) + sqrt(A_slice(1:end-1).*A_slice(2:end))));

%% 6. Method 5 — Delaunay Triangulated Irregular Network (TIN Surface)
tri = delaunay(T.Local_X_m, T.Local_Y_m);
xA = T.Local_X_m(tri(:,1)); yA = T.Local_Y_m(tri(:,1));
xB = T.Local_X_m(tri(:,2)); yB = T.Local_Y_m(tri(:,2));
xC = T.Local_X_m(tri(:,3)); yC = T.Local_Y_m(tri(:,3));
tri_area = 0.5 * abs((xB - xA).*(yC - yA) - (xC - xA).*(yB - yA));
tri_depth = (T.Cut_Depth_m(tri(:,1)) + T.Cut_Depth_m(tri(:,2)) + T.Cut_Depth_m(tri(:,3))) / 3;
V_tin = sum(tri_area .* tri_depth);

%% 7. Print Comparative Summary Table
methods = {'End-Area (Trapezoidal)'; 'Prismoidal (Simpson 1/3)'; ...
           'Borrow-Pit Grid (1-2-3-4)'; 'Contour Slicing (1.0m CI)'; 'Delaunay TIN Surface'};
volumes = [V_end_area; V_prismoidal; V_grid; V_contour; V_tin];
Results = table(methods, volumes, 'VariableNames', {'Method', 'Volume_m3'});
disp(Results);

%% 8. 3D Excavation Surface Plot
figure('Color', 'w', 'Position', [100 100 900 520]);
surf(XX, YY, Z_post, D, 'EdgeAlpha', 0.25);
colormap(turbo); cb = colorbar; cb.Label.String = 'Cut Depth (m)';
xlabel('Easting (m E — UTM Zone 32N)');
ylabel('Northing (m N — UTM Zone 32N)');
zlabel('Post-Excavation Elevation (m)');
title('Abuja Central Area Deep Excavation DEM — 3D Survey Surface');
view(-38, 30); grid on;
