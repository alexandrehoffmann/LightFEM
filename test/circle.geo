lc=0.5; //characteristic mesh size (optional make smaller to refine mesh)

// Place points
Point(1) = { 0, 0, 0, lc};
Point(2) = { 1, 0, 0, lc};
Point(3) = { 0, 1, 0, lc};
Point(4) = {-1, 0, 0, lc};
Point(5) = { 0,-1, 0, lc};

// Create circles from points

// first, center, end
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// create surface
Line Loop(1) = {2, 3, 4, 1};
Plane Surface(1) = {1};

// Showing off another gmsh feature:
// Tell gmsh to use mostly quads rather than tris for surface #2
Recombine Surface(1);

Physical Curve("dOmega") = {2, 1, 4, 3};
Physical Surface("Omega") = {1};
