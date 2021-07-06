//SetFactory("OpenCASCADE");

//Mesh.Algorithm = 6;

//R1 = 90;
//R2 = 63;
//Rn = 0.2*R2;
//h1 = Sqrt(R1*R1 - Rn*Rn);
//h2 = Sqrt(R2*R2 - Rn*Rn);

ratio_lc = 0.1; //0.1 - coarse computations, 0.02 - fine computations
lc = ratio_lc * R1;
lcr = 0.5*lc;

O1 = 0.0;

Point(1) = {O1, 0, 0, lc};
Point(2) = {O1+h1, 0, 0, lcr};
Point(3) = {O1+h1,-Rn, 0, lcr};
Point(4) = {O1+h1, Rn, 0, lcr};
Point(5) = {O1, R1, 0, lc};
Point(6) = {O1-R1, 0, 0, lc};
Point(7) = {O1,-R1, 0, lc};

// Second sphere
O2 = O1 + h1 + h2;
//O2 = O1 + offset;
Point(8) = {O2, 0, 0, lc};
Point(9) = {O2,-R2, 0, lc};
Point(10) = {O2+R2, 0, 0, lc};
Point(11) = {O2, R2, 0, lc};

// Left sphere
Circle(2) = {4,1,5};
Circle(3) = {5,1,6};
Circle(4) = {6,1,7};
Circle(5) = {7,1,3};

Circle(6) = {3,8,9};
Circle(7) = {9,8,10};
Circle(8) = {10,8,11};
Circle(9) = {11,8,4};

// Surfaces
Curve Loop(1) = {2,3,4,5,6,7,8,9};
Plane Surface(1) = {1};

Point{1} In Surface {1};
Point{8} In Surface {1};

Printf("x1 = %g", O1);
Printf("x2 = %g", O2);
