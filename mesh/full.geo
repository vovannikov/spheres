ls = 0.1;

X1 = 1.14;
X2 = 0.85;
X3 = 1.43;

Y1 = 1.35;
Y2 = 0.5;
Y3 = 4.35;

Point(1) = {X1+X2+X3, 0, 0, ls};
Point(2) = {X1+X2, Y3, 0, ls};
Point(3) = {X1+X2, Y2, 0, ls};
Point(4) = {X1, Y2, 0, ls};
Point(5) = {X1, Y1, 0, ls};
Point(6) = {-X1, Y1, 0, ls};
Point(7) = {-X1, Y2, 0, ls};
Point(8) = {-X1-X2, Y2, 0, ls};
Point(9) = {-X1-X2, Y3, 0, ls};
Point(10) = {-X1-X2-X3, 0, 0, ls};
Point(11) = {-X1-X2, -Y3, 0, ls};
Point(12) = {-X1-X2, -Y2, 0, ls};
Point(13) = {-X1, -Y2, 0, ls};
Point(14) = {-X1, -Y1, 0, ls};
Point(15) = {X1, -Y1, 0, ls};
Point(16) = {X1, -Y2, 0, ls};
Point(17) = {X1+X2, -Y2, 0, ls};
Point(18) = {X1+X2, -Y3, 0, ls};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};

Plane Surface(1) = {1};
