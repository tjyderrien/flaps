load 'meshElements.dat'

M=151
N=151

x=meshElements(:,3);
y=meshElements(:,4);

x2=reshape(x,[M,N]);
y2=reshape(y,[M,N]);

%normal importation
Nnx=meshElements(:,5);
Nny=meshElements(:,6);
Nsx=meshElements(:,7);
Nsy=meshElements(:,8);

Nex=meshElements(:,9);
Ney=meshElements(:,10);
Nwx=meshElements(:,11);
Nwy=meshElements(:,12);

Nnx2=reshape(Nnx,[M,N]);
Nny2=reshape(Nny,[M,N]);
Nsx2=reshape(Nsx,[M,N]);
Nsy2=reshape(Nsy,[M,N]);
Nex2=reshape(Nex,[M,N]);
Ney2=reshape(Ney,[M,N]);
Nwx2=reshape(Nwx,[M,N]);
Nwy2=reshape(Nwy,[M,N]);

%area import
AreaN=meshElements(:,13);
AreaS=meshElements(:,14);
AreaE=meshElements(:,15);
AreaW=meshElements(:,16);

AreaN2=reshape(AreaN,[M,N]);
AreaS2=reshape(AreaS,[M,N]);
AreaE2=reshape(AreaE,[M,N]);
AreaW2=reshape(AreaW,[M,N]);

% correction for corners
AreaS2(N,1)=AreaS2(N,1)
AreaE2(N,1)=AreaE2(N,1)
AreaN2(1,1)=AreaN2(1,1)
AreaE2(1,1)=AreaE2(1,1)

%volume import
Vol=meshElements(:,17);
Vol2=reshape(Vol,[M,N]);

%tangent import
TangentNx=meshElements(:,18);
TangentNy=meshElements(:,19);
TangentSx=meshElements(:,20);
TangentSy=meshElements(:,21);
TangentEx=meshElements(:,22);
TangentEy=meshElements(:,23);
TangentWx=meshElements(:,24);
TangentWy=meshElements(:,25);

TangentNx2=reshape(TangentNx,[M,N]);
TangentNy2=reshape(TangentNy,[M,N]);
TangentSx2=reshape(TangentSx,[M,N]);
TangentSy2=reshape(TangentSy,[M,N]);
TangentEx2=reshape(TangentEx,[M,N]);
TangentEy2=reshape(TangentEy,[M,N]);
TangentWx2=reshape(TangentWx,[M,N]);
TangentWy2=reshape(TangentWy,[M,N]);

%test with no tangential transport
% TangentNx2=zeros(M,N); TangentNy2=zeros(M,N)
% TangentSx2=zeros(M,N); TangentSy2=zeros(M,N)
% TangentEx2=zeros(M,N); TangentEy2=zeros(M,N)
% TangentWx2=zeros(M,N); TangentWy2=zeros(M,N)

% irregular mesh directions
CurviNx=meshElements(:,26); CurviNy=meshElements(:,27);
CurviSx=meshElements(:,28); CurviSy=meshElements(:,29);
CurviEx=meshElements(:,30); CurviEy=meshElements(:,31);
CurviWx=meshElements(:,32); CurviWy=meshElements(:,33);

CurviNx2=reshape(CurviNx, [M,N]); 
CurviNy2=reshape(CurviNy, [M,N]);
CurviSx2=reshape(CurviSx, [M,N]);
CurviSy2=reshape(CurviSy, [M,N]);
CurviEx2=reshape(CurviEx, [M,N]);
CurviEy2=reshape(CurviEy, [M,N]);
CurviWx2=reshape(CurviWx, [M,N]);
CurviWy2=reshape(CurviWy, [M,N]);

% distance vectors
DistN=meshElements(:,34);
DistS=meshElements(:,35);
DistE=meshElements(:,36);
DistW=meshElements(:,37);

DistN2=reshape(DistN, [M,N]);
DistS2=reshape(DistS, [M,N]); 
DistE2=reshape(DistE, [M,N]);
DistW2=reshape(DistW, [M,N]); 

% distance entre duals
DistDualN=meshElements(:,38); 
DistDualS=meshElements(:,39);
DistDualE=meshElements(:,40);
DistDualW=meshElements(:,41); 

DistDualN2=reshape(DistDualN, [M,N]);
DistDualS2=reshape(DistDualS, [M,N]); 
DistDualE2=reshape(DistDualE, [M,N]);
DistDualW2=reshape(DistDualW, [M,N]); 

% others

Ex=-1e9
Ey=0e0

ExE=(Ex.*TangentEx2+Ey.*TangentEy2).*TangentEx2;
EyE=(Ex.*TangentEx2+Ey.*TangentEy2).*TangentEy2;
ExN=(Ex.*TangentNx2+Ey.*TangentNy2).*TangentNx2;
EyN=(Ex.*TangentNx2+Ey.*TangentNy2).*TangentNy2;
ExS=(Ex.*TangentSx2+Ey.*TangentSy2).*TangentSx2;
EyS=(Ex.*TangentSx2+Ey.*TangentSy2).*TangentSy2;
ExW=(Ex.*TangentWx2+Ey.*TangentWy2).*TangentWx2;
EyW=(Ex.*TangentWx2+Ey.*TangentWy2).*TangentWy2;

% figure(3)

% quiver(x2(:,1:5),y2(:,1:5),TangentNx2(:,1:5),TangentNy2(:,1:5))
% title 'Tangents'

% quiver(x2,y2,(ExN+ExS+ExE+ExW)/4,(EyN+EyS+EyE+EyW)/4)
% title 'Projected field'

% figure(4)

FluxNbnd=ExN.*TangentNx2+EyN.*TangentNy2;
FluxSbnd=ExS.*TangentSx2+EyS.*TangentSy2;
FluxEbnd=ExE.*TangentEx2+EyE.*TangentEy2;
FluxWbnd=ExW.*TangentWx2+EyW.*TangentWy2;

% FluxNbnd=ExN.*TangentEx2+EyN.*TangentEy2;
% FluxSbnd=ExS.*Nsx2+EyS.*Nsy2;
% FluxEbnd=ExE.*Nex2+EyE.*Ney2;
% FluxWbnd=ExW.*Nwx2+EyW.*Nwy2;

% surf(x2,y2,FluxNbnd+FluxSbnd+FluxEbnd+FluxWbnd)

% flux corrige
FluxN=ExN.*Nnx2+EyN.*Nny2;
FluxS=ExS.*Nsx2+EyS.*Nsy2;
FluxE=ExE.*Nex2+EyE.*Ney2;
FluxW=ExW.*Nwx2+EyW.*Nwy2;


% flux 
% FluxW=(Ex.*TangentWx2+Ey.*TangentWy2).*(TangentWx2.*Nwx2+TangentWy2.*Nwy2) + (Ex.*Nwx2+Ey.*Nwy2).*(Nwx2.^2+Nwy2.^2);
% FluxS=(Ex.*TangentSx2+Ey.*TangentSy2).*(TangentSx2.*Nsx2+TangentSy2.*Nsy2) + (Ex.*Nsx2+Ey.*Nsy2).*(Nsx2.^2+Nsy2.^2);
% FluxE=(Ex.*TangentEx2+Ey.*TangentEy2).*(TangentEx2.*Nex2+TangentEy2.*Ney2) + (Ex.*Nex2+Ey.*Ney2).*(Nex2.^2+Ney2.^2);
% FluxN=(Ex.*TangentNx2+Ey.*TangentNy2).*(TangentNx2.*Nnx2+TangentNy2.*Nny2) + (Ex.*Nnx2+Ey.*Nny2).*(Nnx2.^2+Nny2.^2);

% corrected flux
% FluxW=(Ex.*TangentWx2+Ey.*TangentWy2).*(TangentWx2.^2+TangentWy2.^2) + (Ex.*Nwx2+Ey.*Nwy2).*(Nwx2.^2+Nwy2.^2);
% FluxS=(Ex.*TangentSx2+Ey.*TangentSy2).*(TangentSx2.^2+TangentSy2.^2) + (Ex.*Nsx2+Ey.*Nsy2).*(Nsx2.^2+Nsy2.^2);
% FluxE=(Ex.*TangentEx2+Ey.*TangentEy2).*(TangentEx2.^2+TangentEy2.^2) + (Ex.*Nex2+Ey.*Ney2).*(Nex2.^2+Ney2.^2);
% FluxN=(Ex.*TangentNx2+Ey.*TangentNy2).*(TangentNx2.^2+TangentNy2.^2) + (Ex.*Nnx2+Ey.*Nny2).*(Nnx2.^2+Nny2.^2);

% normal diffusion flux

DenomN=Nnx2.*CurviNx2+Nny2.*CurviNy2;
DenomS=Nsx2.*CurviSx2+Nsy2.*CurviSy2;
DenomE=Nex2.*CurviEx2+Ney2.*CurviEy2;
DenomW=Nwx2.*CurviWx2+Nwy2.*CurviWy2;

DiffNormN=(Nnx2.*Nnx2+Nny2.*Nny2)./(Nnx2.*CurviNx2+Nny2.*CurviNy2);
DiffNormS=(Nsx2.*Nsx2+Nsy2.*Nsy2)./(Nsx2.*CurviSx2+Nsy2.*CurviSy2);
DiffNormE=(Nex2.*Nex2+Ney2.*Ney2)./(Nex2.*CurviEx2+Ney2.*CurviEy2);
DiffNormW=(Nwx2.*Nwx2+Nwy2.*Nwy2)./(Nwx2.*CurviWx2+Nwy2.*CurviWy2);

DiffCrossN=(CurviNx2.*TangentNx2+CurviNy2.*TangentNy2)./(Nnx2.*CurviNx2+Nny2.*CurviNy2);
DiffCrossS=(CurviSx2.*TangentSx2+CurviSy2.*TangentSy2)./(Nsx2.*CurviSx2+Nsy2.*CurviSy2);
DiffCrossE=(CurviEx2.*TangentEx2+CurviEy2.*TangentEy2)./(Nex2.*CurviEx2+Ney2.*CurviEy2);
DiffCrossW=(CurviWx2.*TangentWx2+CurviWy2.*TangentWy2)./(Nwx2.*CurviWx2+Nwy2.*CurviWy2);
