// Length caracteristic
lc = 0.2;
Height=1;
LengthTip=0.2;
LengthBody=2.0;
Length=3.0; 
//Mesh dimensions
NumM=151;
NumN=151;
BumpValueX=0.0001;
BumpValueY=1;

// Summits

apex=newp; Point(apex) = {0,0,0,lc};
NorthControl=newp; Point(NorthControl) = {LengthBody,0.5,0,lc};
SouthControl=newp; Point(SouthControl) = {LengthBody,-0.5,0,lc};

NE=newp; Point(NE) = {Length,Height/2,0,lc};
SE=newp; Point(SE) = {Length,-Height/2,0,lc};

NW=newp; Point(NW) = {LengthTip,0.3,0,lc};
SW=newp; Point(SW) = {LengthTip,-0.3,0,lc};


// Defining contours
North=newl; Spline(North) = {NE, NorthControl, NW};
West=newl; Spline(West) = {NW, apex, SW};
South=newl; Spline(South) = {SW, SouthControl, SE}; 
East=newl; Line(East) = {SE, NE};

//Generate the oriented contour
Contour=newl; Line Loop(Contour) = {North,West,South,East};

// Defining surface from oriented contour
Needle=newl; 
Plane Surface(Needle) = {Contour};

//CONFIG A STRUCTURED MESH
NorthGrid=newl; SouthGrid=NorthGrid+1; Transfinite Line{North,South} = NumM Using Bump BumpValueX; 
WestGrid=newl; EastGrid=WestGrid+1; Transfinite Line{West,East} = NumN Using Bump BumpValueY; 

Transfinite Surface{Needle} = {NE, NW, SW, SE};

Recombine Surface{Needle}; 
Mesh.Smoothing=1000; 

Physical Surface(Needle) = 1;

//SAVING MESH
SetOrder 2;
Save Sprintf("mesh%03g.msh", 0);

