lbathy=0.07; //characteristic mesh size (optional make smaller to refine mesh)
lshallow=0.14	; //characteristic mesh size (optional make smaller to refine mesh)
ldeep=1; //characteristic mesh size (optional make smaller to refine mesh)

zbathy = -0.07;
zshallow  = -2.0;
zdeep = -5.0;

ylim = 15.05;

// Place points

// bathy
Point(1) = {   0,               0, 0};
Point(2) = {ylim,               0, 0};
Point(3) = {ylim,          zbathy, 0, lshallow};
Point(8) = {   0,          zbathy, 0, lshallow};
// shallow 
Point(4) = {ylim,        zshallow, 0, lshallow};
Point(7) = {   0,        zshallow, 0, lshallow};
// deep
Point(5) = {ylim,           zdeep, 0, ldeep};
Point(6) = {   0,           zdeep, 0, ldeep};

// bathy
Line(1) = {8, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 8};

// shallow 
Line(5) = {8, 7};
Line(6) = {7, 4};
Line(7) = {4, 3};

// deep
Line(8) = {7, 6};
Line(9) = {6, 5};
Line(10) = {5, 4};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

Curve Loop(2) = {5, 6, 7, 4};
Plane Surface(2) = {2};

Curve Loop(3) = {8, 9, 10, -6};
Plane Surface(3) = {3};

Physical Surface("bathy") = {1};
Physical Surface("shallow") = {2};
Physical Surface("deep") = {3};

Transfinite Curve {1, 3} = 1 Using Progression 1;
Transfinite Curve {2, 4} = 215 Using Progression 1;

Physical Curve("neumann") = {5, 8, 9, 10, 7, 3, 1};
Physical Curve("dirichlet") = {2};

Recombine Surface(1);
Recombine Surface(2);
Recombine Surface(3);


