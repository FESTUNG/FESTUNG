cx = 0.0;  cy = 0.0;  r = 0.5;
Point(0) = {cx, cy,   0.0, 1.0};  Point(1) = {cx-r, cy, 0.0, 1.0};  Point(2) = {cx, cy+r, 0.0, 1.0};
Point(3) = {cx+r, cy, 0.0, 1.0};  Point(4) = {cx, cy-r, 0.0, 1.0};
Circle(5) = {2, 0, 1};  Circle(6) = {3, 0, 2};  Circle(7) = {4, 0, 3};  Circle(8) = {1, 0, 4};
Line Loop(9) = {5, 6, 7, 8}; // boundary edge IDs will be 5 to 8
Plane Surface(10) = {9};
