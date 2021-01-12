ls = 0.05;

X1 = 1.14;
X2 = 0.85;
X3 = 1.43;

Y1 = 1.35;
Y2 = 0.5;
Y3 = 4.35;

Point(1) = {0, 0, 0, ls};
Point(2) = {X1+X2+X3, 0, 0, ls};
Point(3) = {X1+X2, Y3, 0, ls};
Point(4) = {X1+X2, Y2, 0, ls};
Point(5) = {X1, Y2, 0, ls};
Point(6) = {X1, Y1, 0, ls};
Point(7) = {0, Y1, 0, ls};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};

Plane Surface(1) = {1};
