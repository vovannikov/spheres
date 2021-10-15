//SetFactory("OpenCASCADE");

// This is required to generate better mesh for torsion case
Mesh.Algorithm3D = 10;

//R1 = 90;
//R2 = 63;
//Rn = 0.2*R2;
//h1 = Sqrt(R1*R1 - Rn*Rn);
//h2 = Sqrt(R2*R2 - Rn*Rn);
//ratio_lc = 0.15; // 2.0 - super coarse for debug //0.1 - coarse computations, 0.02 - fine computations

ratio_lc = lc;
lc1 = ratio_lc * R1;
lc2 = ratio_lc * R2;

// Here we assume that R2 <= R1
//lcr = 0.5*Min(lc1, lc2);
lcr = 0.5*lc2;

O1 = 0.0;
Point(1) = {O1, 0, 0, lc1};

Point(2) = {O1+h1, 0, 0, lcr};
Point(3) = {O1+h1,-Rn, 0, lcr};
Point(4) = {O1+h1, Rn, 0, lcr};
Point(5) = {O1+h1, 0,-Rn, lcr};
Point(6) = {O1+h1, 0, Rn, lcr};

Point(7) = {O1, R1, 0, lc1};
Point(8) = {O1-R1, 0, 0, lc1};
Point(9) = {O1,-R1, 0, lc1};
Point(10) = {O1, 0,-R1, lc1};
Point(11) = {O1, 0, R1, lc1};

// Neck
Circle(1) = {3,2,6};
Circle(2) = {6,2,4};
Circle(3) = {4,2,5};
Circle(4) = {5,2,3};

// Left sphere
Circle(5) = {4,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,3};
Circle(9) = {6,1,11};
Circle(10) = {11,1,8};
Circle(11) = {8,1,10};
Circle(12) = {10,1,5};
Circle(13) = {7,1,11};
Circle(14) = {11,1,9};
Circle(15) = {9,1,10};
Circle(16) = {10,1,7};

// Surfaces
Curve Loop(17) = {14,8,1,9};
Surface(18) = {17};
Curve Loop(19) = {13,-9,2,5};
Surface(20) = {19};
Curve Loop(21) = {6,-10,-13};
Surface(22) = {21};
Curve Loop(23) = {7,-14,10};
Surface(24) = {23};
Curve Loop(25) = {15,12,4,-8};
Surface(26) = {25};
Curve Loop(27) = {-16,12,-3,5};
Surface(28) = {27};
Curve Loop(29) = {6,11,16};
Surface(30) = {29};
Curve Loop(31) = {11,-15,-7};
Surface(32) = {31};

// Neck part
Curve Loop(33) = {-1,-2,-3,-4};
Surface(34) = {33};

// Middle of the left sphere
Curve Loop(35) = {13,14,15,16};
Surface(36) = {35};

Surface Loop(71) = {32,30,24,22,36};
Volume(72) = {71};

Surface Loop(73) = {36,28,26,20,18,34};
Volume(74) = {73};

Point{1} In Volume {72};
Point{1} In Volume {74};

// Second sphere
O2 = O1 + h1 + h2;
Point(12) = {O2, 0, 0, lc2};
Point(13) = {O2+R2, 0, 0, lc2};
Point(14) = {O2,-R2, 0, lc2};
Point(15) = {O2, R2, 0, lc2};
Point(16) = {O2, 0,-R2, lc2};
Point(17) = {O2, 0, R2, lc2};

// Right sphere
Circle(37) = {3,12,14};
Circle(38) = {14,12,13};
Circle(39) = {13,12,15};
Circle(40) = {15,12,4};
Circle(41) = {5,12,16};
Circle(42) = {16,12,13};
Circle(43) = {13,12,17};
Circle(44) = {17,12,6};
Circle(45) = {17,12,14};
Circle(46) = {14,12,16};
Circle(47) = {16,12,15};
Circle(48) = {15,12,17};

// Surfaces
Curve Loop(49) = {37,-45,44,-1};
Surface(50) = {49};
Curve Loop(51) = {-44,-48,40,-2};
Surface(52) = {51};
Curve Loop(53) = {38,43,45};
Surface(54) = {53};
Curve Loop(55) = {39,48,-43};
Surface(56) = {55};
Curve Loop(57) = {41,-46,-37,-4};
Surface(58) = {57};
Curve Loop(59) = {40,3,41,47};
Surface(60) = {59};
Curve Loop(61) = {46,42,-38};
Surface(62) = {61};
Curve Loop(63) = {-47,42,39};
Surface(64) = {63};

// Middle of the right sphere
Curve Loop(65) = {45,46,47,48};
Surface(66) = {65};

Surface Loop(75) = {64,62,56,54,66};
Volume(76) = {75};

Surface Loop(77) = {66,60,58,52,50,34};
Volume(78) = {77};

Point{12} In Volume {76};
Point{12} In Volume {78};

// Physical volumes
Physical Volume(81) = {72};
Physical Volume(82) = {74};
Physical Volume(83) = {76};
Physical Volume(84) = {78};

Printf("x1 = %g", O1);
Printf("x2 = %g", O2);