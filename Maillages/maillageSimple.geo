a = 0.0;
b=20.;


Point(1) = {0,0,0,1.0};
Point(2) = {b,0,0,1.0};
Point(3) = {b,b,0,1.0};
Point(4) = {0,b,0,1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};

Plane Surface(6) = {5};

Physical Surface("OMEGA") = {6};
Transfinite Line {1, 2, 3, 4} = 81 Using Progression 1;
Transfinite Surface "*";
