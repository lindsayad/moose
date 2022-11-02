SetFactory ("OpenCASCADE");

r1 = 0.425;
r2 = 0.475;
p = 1.14;
lc = 1;
height = 90;


Point(1) = {0, 0, 0, 1};
Point(2) = {r1*Sqrt(1/2), r1*Sqrt(1/2), 0, lc};
Point(3) = {-r1*Sqrt(0.5), r1*Sqrt(0.5), 0, lc};
Point(4) = {-r1*Sqrt(0.5), -r1*Sqrt(0.5), 0, lc};
Point(5) = {r1*Sqrt(0.5), -r1*Sqrt(0.5), 0, lc};
Point(6) = {r2*Sqrt(0.5), r2*Sqrt(0.5), 0, lc};
Point(7) = {-r2*Sqrt(0.5), r2*Sqrt(0.5), 0, lc};
Point(8) = {-r2*Sqrt(0.5), -r2*Sqrt(0.5), 0, lc};
Point(9) = {r2*Sqrt(0.5), -r2*Sqrt(0.5), 0, lc};
Point(10) = {p/2, p/2, 0, lc};
Point(11) = {-p/2, p/2, 0, lc};
Point(12) = {-p/2, -p/2, 0, lc};
Point(13) = {p/2, -p/2, 0, lc};

//inner boundary
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
//mid boundary
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
//outer boundary
Line(9) = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 13};
Line(12) = {13, 10};
// inner auxiliary line
Line(13) = {1, 2};
Line(14) = {1, 3};
Line(15) = {1, 4};
Line(16) = {1, 5};
//mid auxiliary line
Line(17) = {2, 6};
Line(18) = {3, 7};
Line(19) = {4, 8};
Line(20) = {5, 9};
//outer auxiliary line
Line(21) = {6, 10};
Line(22) = {7, 11};
Line(23) = {8, 12};
Line(24) = {9, 13};
//inner region
Curve Loop(1) = {13, 1, -14};
Plane Surface(1) = {1};
Curve Loop(2) = {14, 2, -15};
Plane Surface(2) = {2};
Curve Loop(3) = {15, 3, -16};
Plane Surface(3) = {3};
Curve Loop(4) = {16, 4, -13};
Plane Surface(4) = {4};
//mid region
Curve Loop(5) = {17, 5, -18, -1};
Plane Surface(5) = {5};
Curve Loop(6) = {18, 6, -19, -2};
Plane Surface(6) = {6};
Curve Loop(7) = {19, 7, -20, -3};
Plane Surface(7) = {7};
Curve Loop(8) = {20, 8, -17, -4};
Plane Surface(8) = {8};
//outer region
Curve Loop(9) = {21, 9, -22, -5};
Plane Surface(9) = {9};
Curve Loop(10) = {22, 10, -23, -6};
Plane Surface(10) = {10};
Curve Loop(11) = {23, 11, -24, -7};
Plane Surface(11) = {11};
Curve Loop(12) = {24, 12, -21, -8};
Plane Surface(12) = {12};

Recombine Surface {1:12};

Transfinite Curve {1:12} = 10 Using Progression 1; // Circumferential layers of 1/4 circle
Transfinite Curve {13:16} = 5 Using Progression 1;// Radial layers of circles
Transfinite Curve {17:20} = 2 Using Progression 1;// Radial layers of mid region
Transfinite Curve {21:24} = 5 Using Progression 1;// Radial layers of outer region

Transfinite Surface {1:12};

Extrude {0, 0, height} {
  Surface{1:12};
  Layers{40};
  Recombine;
}

Physical Point("pin1") = {1};
Physical Point("pin2") = {14};
Physical Point("pin3") = {23};

Physical Surface("inner_bottom") = {1, 2, 3, 4};
Physical Surface("inner_top") = {16, 19, 22, 24};
Physical Surface("inner_wall") = {14, 17, 20, 23};
Physical Surface("mid_bottom") = {5, 6, 7, 8};
Physical Surface("mid_top") = {28, 31, 34, 36};
Physical Surface("mid_wall") = {26, 29, 32, 35};
Physical Surface("outer_bottom") = {9, 10, 11, 12};
Physical Surface("outer_top") = {40, 43, 46, 48};
Physical Surface("outer_wall") = {41, 47, 44, 38};

Physical Volume("inner") = {1, 2, 3, 4};
Physical Volume("mid") = {5, 6, 7, 8};
Physical Volume("outer") = {9, 10, 11, 12};
