//T.J-Y. Derrien, Hilase Centre, Praha
//Thu Apr 9 2015

// /!\ Careful: don't use function names for new variables. But, case is sensitive.

//Contour: length caracteristic
//length unit
lc = 1; 

//cone height
height=1.0;
//cone length
length=5.0;

LengthTip=0.05*length; //X coord of the NW and SW points of the tip
HeightTip=0.05*height; //Y coord of the NW and SW points of the tip

LengthControl1=0.7*length; //X position of the control points for N and S lines
HeightControl1=0.8*height/2; //Y position of the control points for N and S lines

LengthControl2=0.1*length; //X position of the control points for N and S lines
HeightControl2=0.1*height/2; //Y position of the control points for N and S lines

//Mesh dimensions
NumM=401;
NumN=301;

// Mesh refinement on contours
BumpValueX=1.0;
BumpValueY=1.0;
// Mesh.CharacteristicLengthFromPoints=1;
// Mesh.CharacteristicLengthFromCurvature=0;
//GmshSetOption("Mesh", "Algorithm", 5);
//Smoothing and refinement

// SmoothIterations=5;

// CORNERS definition

apex=newp; Point(apex) = {0,0,0, lc};
NorthControl1=newp; Point(NorthControl1) = {LengthControl1,HeightControl1,0,lc};
SouthControl1=newp; Point(SouthControl1) = {LengthControl1,-HeightControl1,0,lc};

NW=newp; Point(NW) = {LengthTip,HeightTip,0,lc};
NE=newp; Point(NE) = {length, height/2, 0, lc};
SE=newp; Point(SE) = {length, -height/2, 0, lc};
SW=newp; Point(SW) = {LengthTip,-HeightTip,0,lc};

// NorthControl2=newp; Point(NorthControl2) = {LengthControl2,HeightControl2,0,lc};
// SouthControl2=newp; Point(SouthControl2) = {LengthControl2,-HeightControl2,0,lc};

// Defining contours
North=newl; Spline(North) = {NW, NorthControl1, NE};
East=newl; Line(East) = {NE, SE};
South=newl; Spline(South) = {SE, SouthControl1, SW};
West=newl; Spline(West) = {SW, apex, NW};

//Delete {NorthControl1; SouthControl1; apex;}; //not effective if linked with elements


//Generate the oriented contour
Contour=newl; Line Loop(Contour) = {North,West,South,East};


// Defining surface from oriented contour
Needle=newl;
Plane Surface(Needle) = {Contour};

//CONFIG A STRUCTURED MESH
NorthGrid=newl; SouthGrid=NorthGrid+1; Transfinite Line{North,South} = NumM Using Bump BumpValueY;
WestGrid=newl; EastGrid=WestGrid+1; Transfinite Line{West,East} = NumN Using Bump BumpValueX;

Transfinite Surface{Needle} = {NW, NE, SE, SW};

Recombine Surface{Needle};
//Mesh.Smoothing=SmoothIterations;

Physical Surface(Needle) = 1;

//SAVING MESH
//SetOrder 3;
Mesh 2;
// AdaptMesh; 
// Mesh.Format=40;
Mesh.LineNumbers=0;
Mesh.PointNumbers=0; 
// Print.GeoOnlyPhysicals=1;
Save Sprintf("mesh.msh");

// Save Sprintf("mesh%03g.msh", 0);

