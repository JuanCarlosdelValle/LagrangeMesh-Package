(* ::Package:: *)

BeginPackage["LagrangeMesh`"];


LagMeshEigenvalues::usage = "LagMeshEigenvalues[] calculates the eigenvalues of the Schr\[ODoubleDot]dinger equation with potential V(x) using the Lagrange Mesh Method.";
LagMeshEigenfunctions::usage = "LagMeshEigenfunctions[] calculates the eigenfunctions of the Schr\[ODoubleDot]dinger equation with potential V(x) using the Lagrange Mesh Method.";
LagMeshEigensystem::usage = "LagMeshEigensystem[] calculates the eigenvalues and eigenfunctions of the Schr\[ODoubleDot]dinger equation with potential V(x) using the Lagrange Mesh Method.";
BuildMesh::usage = "BuildMesh[] constructs the mesh points used by the Lagrange Mesh Method.";
AvailableMeshQ::usage = "AvailableMeshQ[] shows the meshes previously constructed.";



Unprotect[{Scaling,Dimension,Mass,PotentialShift,ExpectationValue,Weights,DiscreteFunction,CoefficientsOnly,PrintMesh,PrintDomain}];


Begin["`Private`"];


(* Main Functions *)

SetAttributes[LagMeshEigenvalues,{HoldAll,ReadProtected}];
SetAttributes[LagMeshEigenfunctions,{HoldAll,ReadProtected}];
SetAttributes[LagMeshEigensystem,{HoldAll,ReadProtected}];
SetAttributes[BuildMesh,{HoldAll,ReadProtected}];
SetAttributes[AvailableMeshQ,{HoldAll,ReadProtected}];



(* ::Section:: *)
(*LagMeshEigenvalues[]*)


(* Block for LagMeshEigenvalues *)

Options[LagMeshEigenvalues]=Join[Options[Eigenvalues],
{WorkingPrecision->MachinePrecision,Scaling->1,Mass->1,PotentialShift->0}];
SyntaxInformation[LagMeshEigenvalues]={"ArgumentPattern"->{_,{_,_,_},_,_,OptionsPattern[]},"LocalVariables"->{"Plot",{1,\[Infinity]}}}





LagMeshEigenvalues[potential_,{variable_,a_,b_},NLevels_,meshpoints_,OptionsPattern[]]:=
Module[

{V,x,type,meshList,meshName,T,U,H,n,h,min,t,y,\[Mu],PS,\[Sigma],Out},

$Assumptions={Element[variable,Reals],variable!=0};

\[Mu]=OptionValue[Mass];
n=meshpoints;
h=OptionValue[Scaling];
PS=OptionValue[PotentialShift];
\[Sigma]=1;


(* Begin[Alerts] *)
LagMeshEigenvalues::nnarg0 = "The potential `1` is not a function of `2`.";
LagMeshEigenvalues::nnarg1 = "The desired number of levels `1` is not a positive integer number.";
LagMeshEigenvalues::nnarg2 = "Invalid limit(s) in variable `1`. Limits should be real and of the form {`1`, a, b} with a < b.";
LagMeshEigenvalues::nnarg3 = "Value of option `2` -> `1` is not a positive integer number.";
LagMeshEigenvalues::nnarg4 = "Invalid limit(s) in variable `2`. The argumert `1` is not a real number.";
LagMeshEigenvalues::nnarg5 = "Value of option `2` -> `1` is not a real number.";
LagMeshEigenvalues::nnarg6 = "Value of option `2` -> `1` is not positive number.";
LagMeshEigenvalues::nnarg7 = "Value of option `2` -> `1` is not positive number.";
LagMeshEigenvalues::nnarg8 = "The dimension of the mesh should be larger or equal than `1`.";
LagMeshEigenvalues::nnarg9 = "Scaling is not available.";

(* End[Alerts] *)

(* Begin[Abort] *)
If[MemberQ[potential,variable,{0,-1}]==False,If[NumericQ[potential]==False,Message[LagMeshEigenvalues::nnarg0,potential,variable];Abort[]]];
If[Cases[potential,Except[_?NumericQ|variable],{-1}]=={},Continue,Message[LagMeshEigenvalues::nnarg0,potential,variable];Abort[]];

If[SameQ[Simplify[Im[a]],0]==False,Message[LagMeshEigenvalues::nnarg4,a,variable];Abort[]];
If[SameQ[Simplify[Im[b]],0]==False,Message[LagMeshEigenvalues::nnarg4,b,variable];Abort[]];
If[a>=b,Message[LagMeshEigenvalues::nnarg2, variable];Abort[]];
If[NumberQ[a]==False, If[SameQ[a,-\[Infinity]]==False,Message[LagMeshEigenvalues::nnarg4,a,variable];Abort[]]];
If[NumberQ[b]==False, If[SameQ[b,\[Infinity]]==False,Message[LagMeshEigenvalues::nnarg4,b,variable];Abort[]]];
If[IntegerQ[NLevels]==False\[Or]NLevels<=0,Message[LagMeshEigenvalues::nnarg1, NLevels];Abort[]];
If[NLevels>meshpoints,Message[LagMeshEigenvalues::nnarg8,NLevels];Abort[]];



If[
NumericQ[OptionValue[WorkingPrecision]]==False, 
Message[LagMeshEigenvalues::nnarg3, OptionValue[WorkingPrecision],WorkingPrecision];Abort[]];
If[
OptionValue[WorkingPrecision]!=MachinePrecision,
If[
IntegerQ[OptionValue[WorkingPrecision]]==False\[Or]OptionValue[WorkingPrecision]<1, 
Message[LagMeshEigenvalues::nnarg3, OptionValue[WorkingPrecision],WorkingPrecision];Abort[]]];
(*If[NumericQ[OptionValue[Mass]]==False\[Or]SameQ[Simplify[Im[OptionValue[Mass]]],0]==False\[Or]OptionValue[Mass]<=0,Message[LagMeshEigenvalues::nnarg6, OptionValue[Mass],Mass];Abort[]];*)
If[IntegerQ[NLevels]==False\[Or]NLevels<0,Message[LagMeshEigenvalues::nnarg1, NLevels];Abort[]];
If[NumericQ[OptionValue[PotentialShift]]==False\[Or]SameQ[Simplify[Im[OptionValue[PotentialShift]]],0]==False,Message[LagMeshEigenvalues::nnarg5, OptionValue[PotentialShift],PotentialShift];Abort[]];

If[NumericQ[h]==False\[Or]SameQ[Simplify[Im[h]],0]==False\[Or]h<=0,Message[LagMeshEigenvalues::nnarg6, OptionValue[Scaling],Scaling];Abort[]];

(* End[Abort] *)


If[a==-\[Infinity]&&b==\[Infinity],type="Hermite"];
If[a>-\[Infinity]&&b<\[Infinity],type="Legendre"];
If[a>-\[Infinity]&&b==\[Infinity],type="Laguerre";\[Sigma]=1];
If[a==-\[Infinity]&&b<\[Infinity],type="Laguerre";\[Sigma]=-1];

(* Alert for Legendre Mesh *)
If[type=="Legendre",
If[OptionValue[Scaling]!=1,
Message[LagMeshEigenvalues::nnarg9];Abort[]]
];





SetDirectory[NotebookDirectory[]];
Quiet[dir=Directory[];
dir=ToString[dir]<>"/Meshes/"<>type<>"/MeshPoints/";
SetDirectory[dir]];

meshName=type<>"_"<>ToString[n]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>".dat";

If[FileExistsQ[meshName],meshList=Import[meshName,"List"],
BuildMesh[type,n,WorkingPrecision->OptionValue[WorkingPrecision]]];

SetDirectory[dir];
meshList=Import[meshName,"List"];


(* Beginning  of Hermite Case *)
If[type=="Hermite",
V[x_]:=potential/.variable->x;
Do[x[i]=meshList[[i]],{i,1,Length[meshList]}];
(* Kinetic Elements *)
Do[For[i=1,i<j,i++,Subscript[T, i,j]=(-1)^(i-j) 2/(x[i]-x[j])^2;Subscript[T, j,i]=Subscript[T, i,j]],{j,1,n}];
Clear[i];

Do[Subscript[T, i,i]=2/6 (2n+1-x[i]^2),{i,1,n}]

];

(* End of Hermite Case *)


(* Beginning of Lengedre Case *)

If[type=="Legendre",
V[x_]:=potential/.variable->x;
t[y_]:=1/2 (a+b)+1/2 (b-a)y;
meshList=t[meshList];
Do[x[i]=meshList[[i]],{i,1,Length[meshList]}];

Do[For[i=1,i<j,i++,Subscript[T, i,j]=((-1)^(i+j) (a(x[i]+x[j])+b(x[i]+x[j])-2x[i]x[j]-2 a b))/(\[Sqrt]((x[i]-a)(b-x[i]))\[Sqrt]((x[j]-a)(b-x[j])) (x[i]-x[j])^2);Subscript[T, j,i]=Subscript[T, i,j]],{j,1,n}];

Do[Subscript[T, i,i]=((a-b)^2+n(n+1)(x[i]-a)(b-x[i]))/(3(a-x[i])^2 (b-x[i])^2),{i,1,n}];
Clear[i]
];


(* Beginning of Laguerre Case *)

If[type=="Laguerre",
If[a<\[Infinity],
V[x_]:=potential/.variable->x+a];
If[b<\[Infinity],
V[x_]:=potential/.variable->x+b];
Do[x[i]= meshList[[i]],{i,1,Length[meshList]}];
Do[For[i=1,i<j,i++,Subscript[T, i,j]=(-1)^(i-j) (x[i]+x[j])/((x[i]x[j])^(1/2) (x[i]-x[j])^2);Subscript[T, j,i]=Subscript[T, i,j]],{j,1,n}];
Do[Subscript[T, i,i]=(4+(4n+2)x[i]-x[i]^2)/(12 x[i]^2),{i,1,n}];
Clear[i]

];

min=0;
Quiet[If[Abs[
NMinimize[{V[x],a<=x<=b},x,WorkingPrecision->10][[1]]]<\[Infinity],
min=Rationalize[NMinimize[{V[x],a<=x<=b},x,WorkingPrecision->OptionValue[WorkingPrecision]][[1]],10^-OptionValue[WorkingPrecision]],
min=0]];

(* Matrix Potential Elements *)
Quiet[
Do[Subscript[V,i,i]=V[\[Sigma] h x[i]],{i,1,n}]];
U=DiagonalMatrix[Table[Subscript[V,i,i],{i,1,n}]];


(* Hamiltonian Matrix *)
H=Table[Table[1/(2 \[Mu] h^2) Subscript[T,l,m],{m,1,n}],{l,1,n}]+U;


(* Diagonalization *)
Out=Sort[Eigenvalues[H-(min-PS)IdentityMatrix[n],-NLevels,Method->OptionValue[Method],Cubics->OptionValue[Cubics],Quartics->OptionValue[Quartics]]+min-PS];
If[NLevels==1,Out=Out[[1]]];
Out









]


(* ::Section:: *)
(*LagMeshEigenfunctions[]*)


(* Block for LagMeshEigenfunctions *)

Options[LagMeshEigenfunctions]=Join[Options[Eigenvalues],
{WorkingPrecision->MachinePrecision,Scaling->1,Mass->1,PotentialShift->0,ExpectationValue->None,DiscreteFunction->False,CoefficientsOnly->False}];
SyntaxInformation[LagMeshEigenfunctions]={"ArgumentPattern"->{_,{_,_,_},_,_,OptionsPattern[]},"LocalVariables"->{"Plot",{1,\[Infinity]}}}



LagMeshEigenfunctions[potential_,{variable_,a_,b_},NLevels_,meshpoints_,OptionsPattern[]]:=
Module[

{V,x,type,meshList,meshName,weightsList,weightsName,wflag,T,U,H,n,h,min,t,y,\[Mu],f,HermiteN,HH,LegendreN,LaguerreN,vectors,\[Omega],PS,\[Sigma],shift,\[Lambda],fac,EVpotential,V2,EVlst,Out},

$Assumptions=Element[variable,Reals];

\[Mu]=OptionValue[Mass];
n=meshpoints;
h=OptionValue[Scaling];
PS=OptionValue[PotentialShift];
\[Sigma]=1;
shift=0;
EVpotential=OptionValue[ExpectationValue];
V2[x_]:=EVpotential/.variable->x;



(* Begin[Alerts] *)
LagMeshEigenfunctions::nnarg0 = "The potential `1` is not a function of `2`.";
LagMeshEigenfunctions::nnarg1 = "The desired number of levels `1` is not a positive integer number.";
LagMeshEigenfunctions::nnarg2 = "Invalid limit(s) in variable `1`. Limits should be real and of the form {`1`, a, b} with a < b.";
LagMeshEigenfunctions::nnarg3 = "Value of option `2` -> `1` is not a positive integer number.";
LagMeshEigenfunctions::nnarg4 = "Invalid limit(s) in variable `2`. The argumert `1` is not a real number.";
LagMeshEigenfunctions::nnarg5 = "Value of option `2` -> `1` is not a real number.";
LagMeshEigenfunctions::nnarg6 = "Value of option `2` -> `1` is not positive number.";
LagMeshEigenfunctions::nnarg7 = "Value of option `1` -> `2` is not a function of `3`.";
LagMeshEigenfunctions::nnarg8 = "The dimension of the mesh should be larger or equal than `1` .";
LagMeshEigenfunctions::nnarg9 = "Scaling is not available.";
LagMeshEigenfunctions::nnarg10 = "The value of options `2` and `1` cannot be simultaneously True.";
LagMeshEigenfunctions::nnarg11 = "The value of option `1` should be either False or True.";
LagMeshEigenfunctions::nnarg12 = "Value of option ExpectationValue -> `1` is not a function of `2`.";
(* End[Alerts] *)





(* Begin[Abort] *)
If[MemberQ[potential,variable,{0,-1}]==False,If[NumericQ[potential]==False,Message[LagMeshEigenfunctions::nnarg0,potential,variable];Abort[]]];
If[Cases[potential,Except[_?NumericQ|variable],{-1}]=={},Continue,Message[LagMeshEigenfunctions::nnarg0,potential,variable];Abort[]];


If[SameQ[Simplify[Im[a]],0]==False,Message[LagMeshEigenfunctions::nnarg4,a,variable];Abort[]];
If[SameQ[Simplify[Im[b]],0]==False,Message[LagMeshEigenfunctions::nnarg4,b,variable];Abort[]];
If[a>=b,Message[LagMeshEigenfunctions::nnarg2, variable];Abort[]];
If[NumberQ[a]==False, If[SameQ[a,-\[Infinity]]==False,Message[LagMeshEigenfunctions::nnarg4,a,variable];Abort[]]];
If[NumberQ[b]==False, If[SameQ[b,\[Infinity]]==False,Message[LagMeshEigenfunctions::nnarg4,b,variable];Abort[]]];
If[IntegerQ[NLevels]==False\[Or]NLevels<=0,Message[LagMeshEigenfunctions::nnarg1, NLevels];Abort[]];
If[NLevels>meshpoints,Message[LagMeshEigenfunctions::nnarg8,NLevels];Abort[]];
If[
NumericQ[OptionValue[WorkingPrecision]]==False, 
Message[LagMeshEigenfunctions::nnarg3, OptionValue[WorkingPrecision],WorkingPrecision];Abort[]];
If[
OptionValue[WorkingPrecision]!=MachinePrecision,
If[
IntegerQ[OptionValue[WorkingPrecision]]==False\[Or]OptionValue[WorkingPrecision]<1, 
Message[LagMeshEigenfunctions::nnarg3, OptionValue[WorkingPrecision],WorkingPrecision];Abort[]]];



If[NumericQ[OptionValue[Mass]]==False\[Or]SameQ[Simplify[Im[OptionValue[Mass]]],0]==False\[Or]OptionValue[Mass]<=0,Message[LagMeshEigenfunctions::nnarg6, OptionValue[Mass],Mass];Abort[]];
If[IntegerQ[NLevels]==False\[Or]NLevels<0,Message[LagMeshEigenfunctions::nnarg1, NLevels];Abort[]];
If[NumericQ[OptionValue[PotentialShift]]==False\[Or]SameQ[Simplify[Im[OptionValue[PotentialShift]]],0]==False,Message[LagMeshEigenfunctions::nnarg5, OptionValue[PotentialShift],PotentialShift];Abort[]];
If[NumericQ[h]==False\[Or]SameQ[Simplify[Im[h]],0]==False\[Or]h<=0,Message[LagMeshEigenfunctions::nnarg6, OptionValue[Scaling],Scaling];Abort[]];
If[SameQ[OptionValue[ExpectationValue],None]==False,
If[MemberQ[EVpotential,variable,{0,-1}]==False,If[NumericQ[EVpotential]==False,Message[LagMeshEigenfunctions::nnarg7,ExpectationValue,OptionValue[ExpectationValue],variable];Abort[]]],Continue];

If[OptionValue[DiscreteFunction]==True&&OptionValue[CoefficientsOnly]==True, Message[LagMeshEigenfunctions::nnarg10, DiscreteFunction, CoefficientsOnly];Abort[]];
If[MemberQ[{False,True},OptionValue[DiscreteFunction]]==False,Message[LagMeshEigenfunctions::nnarg11, DiscreteFunction];Abort[]];
If[MemberQ[{False,True},OptionValue[CoefficientsOnly]]==False,Message[LagMeshEigenfunctions::nnarg11, CoefficientsOnly];Abort[]];


If[SameQ[EVpotential,None]==False,
If[Cases[EVpotential,Except[_?NumericQ|variable],{-1}]=={},Continue,Message[LagMeshEigenfunctions::nnarg12,EVpotential,variable];Abort[]]];

(* End[Abort] *)


If[a==-\[Infinity]&&b==\[Infinity],type="Hermite"];
If[a>-\[Infinity]&&b<\[Infinity],type="Legendre"];
If[a>-\[Infinity]&&b==\[Infinity],type="Laguerre";\[Sigma]=1;shift=-a];
If[a==-\[Infinity]&&b<\[Infinity],type="Laguerre";\[Sigma]=-1;shift=-b];

(* Alert for Legendre Mesh *)
If[type=="Legendre",
If[OptionValue[Scaling]!=1,
Message[LagMeshEigenfunctions::nnarg9];Abort[]]
];



SetDirectory[NotebookDirectory[]];
dir=Directory[];
Quiet[dir=ToString[dir]<>"/Meshes/"<>type<>"/MeshPoints/";
SetDirectory[dir]];

meshName=type<>"_"<>ToString[n]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>".dat";

If[FileExistsQ[meshName],meshList=Import[meshName,"List"],
BuildMesh[type,n,WorkingPrecision->OptionValue[WorkingPrecision],Weights->OptionValue[DiscreteFunction]]];

SetDirectory[dir];
meshList=Import[meshName,"List"];

ResetDirectory[];
SetDirectory[NotebookDirectory[]];
dir=Directory[];
dir=ToString[dir]<>"/Meshes/"<>type<>"/Weights/";
SetDirectory[dir];
weightsName=type<>"_"<>ToString[n]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>"_Weights.dat";

(* Flag for Weights *)
wflag=0;
If[OptionValue[DiscreteFunction]==True,
If[FileExistsQ[weightsName]==True,wflag=1;
weightsList=Import[weightsName,"List"];
Do[\[Lambda][i]=weightsList[[i]],{i,1,Length[meshList]}]
]];

(* Beginning  of Hermite Case *)

If[type=="Hermite",
V[x_]:=potential/.variable->x;
Do[x[i]=meshList[[i]],{i,1,Length[meshList]}];

(* Kinetic Elements *)
Do[For[i=1,i<j,i++,Subscript[T, i,j]=(-1)^(i-j) 2/(x[i]-x[j])^2;Subscript[T, j,i]=Subscript[T, i,j]],{j,1,n}];
Clear[i];

Do[Subscript[T, i,i]=2/6 (2n+1-x[i]^2),{i,1,n}];

(* EigenFunctions *)
HH=Sqrt[\[Pi]] 2^n n!;
HermiteN[y_]:=2^n \!\(
\*UnderoverscriptBox[\(\[Product]\), \(k = 1\), \(n\)]\((y - x[k])\)\);
f[j_,y_]:=(-1)^(n-j) (2HH)^(-1/2) HermiteN[y]/(y-x[j]);
\[Omega][y_]:=Exp[-(y^2/2)];

If[OptionValue[DiscreteFunction]==True,
(* Weights *)
If[
wflag==0,
Do[\[Lambda][i]=2HH Exp[x[i]^2]/(D[HermiteN[y],y])^2/.y->x[i],{i,1,n}];
weightsList=Table[\[Lambda][i],{i,1,n}];
Export[weightsName,weightsList];
Print[weightsName]];
wflag=1
];



];
(* End of Hermite Case *)


(* Beginning of Lengedre Case *)

If[type=="Legendre",
V[x_]:=potential/.variable->x;
t[y_]:=1/2 (a+b)+1/2 (b-a)y;
meshList=t[meshList];
Do[x[i]=meshList[[i]],{i,1,Length[meshList]}];

Do[For[i=1,i<j,i++,Subscript[T, i,j]=((-1)^(i+j) (a(x[i]+x[j])+b(x[i]+x[j])-2x[i]x[j]-2 a b))/(\[Sqrt]((x[i]-a)(b-x[i]))\[Sqrt]((x[j]-a)(b-x[j])) (x[i]-x[j])^2);Subscript[T, j,i]=Subscript[T, i,j]],{j,1,n}];

Do[Subscript[T, i,i]=((a-b)^2+n(n+1)(x[i]-a)(b-x[i]))/(3(a-x[i])^2 (b-x[i])^2),{i,1,n}];
Clear[i];


(* EigenFunctions *)
LegendreN[y_]:=(1/(b-a))^n (2n)!/(n!)^2 \!\(
\*UnderoverscriptBox[\(\[Product]\), \(k = 1\), \(n\)]\((y - x[k])\)\);
f[j_,y_]:=1/Sqrt[2] (-1)^(j+1) 1/(Sqrt[(x[j]-a)(b-x[j])]) LegendreN[y]/(y-x[j]);
\[Omega][y_]:=((y-a)(b-y))Sqrt[2/(b-a)];
If[OptionValue[DiscreteFunction]==True,
(* Weights *)
If[
wflag==0,
Do[\[Lambda][i]=2/((x[i]-a)(b-x[i])(D[LegendreN[y],y])^2)/.y->x[i],{i,1,n}];
weightsList=Table[\[Lambda][i],{i,1,n}];
Export[weightsName,weightsList];
Print[weightsName];
wflag=1
]
];

];
(* End of Legendre Case *)



(* Beginning of Laguerre Case *)
If[type=="Laguerre",
If[a<\[Infinity],
V[x_]:=potential/.variable->x+a;
V2[x_]:=EVpotential/.variable->x+a];
If[b<\[Infinity],
V[x_]:=potential/.variable->x+b;
V2[x_]:=EVpotential/.variable->x+b];
Do[x[i]=meshList[[i]],{i,1,Length[meshList]}];

Do[For[i=1,i<j,i++,Subscript[T, i,j]=(-1)^(i-j) (x[i]+x[j])/((x[i]x[j])^(1/2) (x[i]-x[j])^2);Subscript[T, j,i]=Subscript[T, i,j]],{j,1,n}];
Do[Subscript[T, i,i]=(4+(4n+2)x[i]-x[i]^2)/(12 x[i]^2),{i,1,n}];

Clear[i];

(* EigenFunctions *)
LaguerreN[y_]:= 1/n! \!\(
\*UnderoverscriptBox[\(\[Product]\), \(k = 1\), \(n\)]\((y - x[k])\)\);
f[j_,y_]:= (-1)^(j+n) x[j]^(-1/2) y LaguerreN[y]/(y-x[j]);
\[Omega][y_]:=Exp[-\[Sigma]((y)/2)];


If[OptionValue[DiscreteFunction]==True,
(* Weights *)
If[
wflag==0,
Do[\[Lambda][i]=Exp[x[i]]/(x[i](D[LaguerreN[y],y])^2)/.y->x[i],{i,1,n}];
weightsList=Table[\[Lambda][i],{i,1,n}];
Export[weightsName,weightsList];
Print[weightsName];
wflag=1]
];

];


(* End of Laguerre Case *)

min=0;
Quiet[If[Abs[
NMinimize[{V[x],a<=x<=b},x,WorkingPrecision->10][[1]]]<\[Infinity],
min=Rationalize[NMinimize[{V[x],a<=x<=b},x,WorkingPrecision->OptionValue[WorkingPrecision]][[1]],10^-OptionValue[WorkingPrecision]],
min=0]];

(* Matrix Potential Elements *)
Quiet[
Do[Subscript[V,i,i]=V[\[Sigma] h x[i]],{i,1,n}]];
U=DiagonalMatrix[Table[Subscript[V,i,i],{i,1,n}]];


(* Hamiltonian Matrix *)
H=Table[Table[1/(2 \[Mu] h^2) Subscript[T,l,m],{m,1,n}],{l,1,n}]+U;



(* Diagonalization *)

vectors=Eigenvectors[H-(min-PS)IdentityMatrix[n],-NLevels,Method->OptionValue[Method],Cubics->OptionValue[Cubics],Quartics->OptionValue[Quartics]];


Do[Do[\!\(\*OverscriptBox[
SubscriptBox[\(c\), \(j\)], \(k\)]\)=vectors[[k]][[j]],{j,1,n}],{k,1,NLevels}];

(* Continuous Wave function *)

If[wflag==0,
If[OptionValue[CoefficientsOnly]==False,
Out=Table[h^(-(1/2)) \[Omega][(variable+shift)/h]SetPrecision[Expand[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(
\*OverscriptBox[
SubscriptBox[\(c\), \(i\)], \(NLevels - k\)] f[i, 
\*FractionBox[\(\[Sigma] \((variable + shift)\)\), \(h\)]]\)\)],OptionValue[WorkingPrecision]],{k,0,NLevels-1}],
Out=Table[Table[\!\(\*OverscriptBox[
SubscriptBox[\(c\), \(i\)], \(NLevels - k\)]\),{i,1,n}],{k,0,NLevels-1}]]
];

(* Discrete Wave function *)
If[wflag==1,
If[type=="Legendre",fac=(-1)^(n+1) Sqrt[2/(b-a)],fac=1];
Out=Table[Table[{\[Sigma] h x[i]-h shift,fac/h^(1/2) \!\(\*OverscriptBox[
SubscriptBox[\(c\), \(i\)], \(NLevels - k\)]\)/Sqrt[\[Lambda][i]]},{i,1,n}],{k,0,NLevels-1}]];
If[NLevels==1,Out=Out[[1]]];

(* Expectation Value *)
If[SameQ[EVpotential,None]==False,
EVlst=Table[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(
\*SuperscriptBox[\((
\*OverscriptBox[
SubscriptBox[\(c\), \(i\)], \(NLevels - k\)])\), \(2\)] V2[\[Sigma]\ x[i]]\)\),{k,0,NLevels-1}];
If[NLevels==1,
Out={Out,EVlst[[1]]},
Out=Transpose[{Out,EVlst}]]
];
Out
]


(* ::Section:: *)
(*LagMeshEigensystem[]*)


(* Block for LagMeshEigensystem *)
Options[LagMeshEigensystem]=Join[Options[Eigenvalues],
{WorkingPrecision->MachinePrecision,Scaling->1,Mass->1,PotentialShift->0,ExpectationValue->None,DiscreteFunction->False,CoefficientsOnly->False}];
SyntaxInformation[LagMeshEigensystem]={"ArgumentPattern"->{_,{_,_,_},_,_,OptionsPattern[]},"LocalVariables"->{"Plot",{1,\[Infinity]}}}



LagMeshEigensystem[potential_,{variable_,a_,b_},NLevels_,meshpoints_,OptionsPattern[]]:=
Module[

{V,x,type,meshList,meshName,weightsList,weightsName,wflag,T,U,H,n,h,min,t,y,\[Mu],f,HermiteN,HH,LegendreN,LaguerreN,vectors,\[Omega],PS,\[Sigma],shift,\[Lambda],fac,EVpotential,V2,EVlst,spectrum,energies,wfunctions,Out},

$Assumptions=Element[variable,Reals];

\[Mu]=OptionValue[Mass];
n=meshpoints;
h=OptionValue[Scaling];
PS=OptionValue[PotentialShift];
\[Sigma]=1;
shift=0;
EVpotential=OptionValue[ExpectationValue];
V2[x_]:=EVpotential/.variable->x;

(* Begin[Alerts] *)
LagMeshEigensystem::nnarg0 = "The potential `1` is not a function of `2`.";
LagMeshEigensystem::nnarg1 = "The desired number of levels `1` is not a positive integer number.";
LagMeshEigensystem::nnarg2 = "Invalid limit(s) in variable `1`. Limits should be real and of the form {`1`, a, b} with a < b.";
LagMeshEigensystem::nnarg3 = "Value of option `2` -> `1` is not a positive integer number.";
LagMeshEigensystem::nnarg4 = "Invalid limit(s) in variable `2`. The argumert `1` is not a real number.";
LagMeshEigensystem::nnarg5 = "Value of option `2` -> `1` is not a real number.";
LagMeshEigensystem::nnarg6 = "Value of option `2` -> `1` is not positive number.";
LagMeshEigensystem::nnarg7 = "Value of option `1` -> `2` is not a function of `3`.";
LagMeshEigensystem::nnarg8 = "The dimension of the mesh should be larger or equal than `1` .";
LagMeshEigensystem::nnarg9 = "Scaling is not available.";
LagMeshEigensystem::nnarg10 = "The value of options `2` and `1` cannot be simultaneously True.";
LagMeshEigensystem::nnarg11 = "The value of option `1` should be either False or True.";
LagMeshEigensystem::nnarg12 = "Value of option ExpectationValue -> `1` is not a function of `2`.";
(* End[Alerts] *)





(* Begin[Abort] *)
If[MemberQ[potential,variable,{0,-1}]==False,If[NumericQ[potential]==False,Message[LagMeshEigensystem::nnarg0,potential,variable];Abort[]]];
If[Cases[potential,Except[_?NumericQ|variable],{-1}]=={},Continue,Message[LagMeshEigensystem::nnarg0,potential,variable];Abort[]];



If[SameQ[Simplify[Im[a]],0]==False,Message[LagMeshEigensystem::nnarg4,a,variable];Abort[]];
If[SameQ[Simplify[Im[b]],0]==False,Message[LagMeshEigensystem::nnarg4,b,variable];Abort[]];
If[a>=b,Message[LagMeshEigensystem::nnarg2, variable];Abort[]];
If[NumberQ[a]==False, If[SameQ[a,-\[Infinity]]==False,Message[LagMeshEigensystem::nnarg4,a,variable];Abort[]]];
If[NumberQ[b]==False, If[SameQ[b,\[Infinity]]==False,Message[LagMeshEigensystem::nnarg4,b,variable];Abort[]]];
If[IntegerQ[NLevels]==False\[Or]NLevels<=0,Message[LagMeshEigensystem::nnarg1, NLevels];Abort[]];
If[NLevels>meshpoints,Message[LagMeshEigensystem::nnarg8,NLevels];Abort[]];
If[
NumericQ[OptionValue[WorkingPrecision]]==False, 
Message[LagMeshEigensystem::nnarg3, OptionValue[WorkingPrecision],WorkingPrecision];Abort[]];
If[
OptionValue[WorkingPrecision]!=MachinePrecision,
If[
IntegerQ[OptionValue[WorkingPrecision]]==False\[Or]OptionValue[WorkingPrecision]<1, 
Message[LagMeshEigensystem::nnarg3, OptionValue[WorkingPrecision],WorkingPrecision];Abort[]]];



If[NumericQ[OptionValue[Mass]]==False\[Or]SameQ[Simplify[Im[OptionValue[Mass]]],0]==False\[Or]OptionValue[Mass]<=0,Message[LagMeshEigensystem::nnarg6, OptionValue[Mass],Mass];Abort[]];
If[IntegerQ[NLevels]==False\[Or]NLevels<0,Message[LagMeshEigensystem::nnarg1, NLevels];Abort[]];
If[NumericQ[OptionValue[PotentialShift]]==False\[Or]SameQ[Simplify[Im[OptionValue[PotentialShift]]],0]==False,Message[LagMeshEigensystem::nnarg5, OptionValue[PotentialShift],PotentialShift];Abort[]];
If[NumericQ[h]==False\[Or]SameQ[Simplify[Im[h]],0]==False\[Or]h<=0,Message[LagMeshEigensystem::nnarg6, OptionValue[Scaling],Scaling];Abort[]];
If[SameQ[OptionValue[ExpectationValue],None]==False,
If[MemberQ[EVpotential,variable,{0,-1}]==False,If[NumericQ[EVpotential]==False,Message[LagMeshEigensystem::nnarg7,ExpectationValue,OptionValue[ExpectationValue],variable];Abort[]]],Continue];

If[OptionValue[DiscreteFunction]==True&&OptionValue[CoefficientsOnly]==True, Message[LagMeshEigensystem::nnarg10, DiscreteFunction, CoefficientsOnly];Abort[]];
If[MemberQ[{False,True},OptionValue[DiscreteFunction]]==False,Message[LagMeshEigensystem::nnarg11, DiscreteFunction];Abort[]];
If[MemberQ[{False,True},OptionValue[CoefficientsOnly]]==False,Message[LagMeshEigensystem::nnarg11, CoefficientsOnly];Abort[]];



If[SameQ[EVpotential,None]==False,
If[Cases[EVpotential,Except[_?NumericQ|variable],{-1}]=={},Continue,Message[LagMeshEigensystem::nnarg12,EVpotential,variable];Abort[]]];
(* End[Abort] *)


If[a==-\[Infinity]&&b==\[Infinity],type="Hermite"];
If[a>-\[Infinity]&&b<\[Infinity],type="Legendre"];
If[a>-\[Infinity]&&b==\[Infinity],type="Laguerre";\[Sigma]=1;shift=-a];
If[a==-\[Infinity]&&b<\[Infinity],type="Laguerre";\[Sigma]=-1;shift=-b];

(* Alert for Legendre Mesh *)
If[type=="Legendre",
If[OptionValue[Scaling]!=1,
Message[LagMeshEigensystem::nnarg9];Abort[]]
];



SetDirectory[NotebookDirectory[]];
dir=Directory[];
Quiet[dir=ToString[dir]<>"/Meshes/"<>type<>"/MeshPoints/";
SetDirectory[dir]];

meshName=type<>"_"<>ToString[n]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>".dat";

If[FileExistsQ[meshName],meshList=Import[meshName,"List"],
BuildMesh[type,n,WorkingPrecision->OptionValue[WorkingPrecision],Weights->OptionValue[DiscreteFunction]]];

SetDirectory[dir];
meshList=Import[meshName,"List"];

ResetDirectory[];
SetDirectory[NotebookDirectory[]];
dir=Directory[];
dir=ToString[dir]<>"/Meshes/"<>type<>"/Weights/";
SetDirectory[dir];
weightsName=type<>"_"<>ToString[n]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>"_Weights.dat";

(* Flag for Weights *)
wflag=0;
If[OptionValue[DiscreteFunction]==True,
If[FileExistsQ[weightsName]==True,wflag=1;
weightsList=Import[weightsName,"List"];
Do[\[Lambda][i]=weightsList[[i]],{i,1,Length[meshList]}]
]];

(* Beginning  of Hermite Case *)

If[type=="Hermite",
V[x_]:=potential/.variable->x;
Do[x[i]=meshList[[i]],{i,1,Length[meshList]}];

(* Kinetic Elements *)
Do[For[i=1,i<j,i++,Subscript[T, i,j]=(-1)^(i-j) 2/(x[i]-x[j])^2;Subscript[T, j,i]=Subscript[T, i,j]],{j,1,n}];
Clear[i];

Do[Subscript[T, i,i]=2/6 (2n+1-x[i]^2),{i,1,n}];

(* EigenFunctions *)
HH=Sqrt[\[Pi]] 2^n n!;
HermiteN[y_]:=2^n \!\(
\*UnderoverscriptBox[\(\[Product]\), \(k = 1\), \(n\)]\((y - x[k])\)\);
f[j_,y_]:=(-1)^(n-j) (2HH)^(-1/2) HermiteN[y]/(y-x[j]);
\[Omega][y_]:=Exp[-(y^2/2)];

If[OptionValue[DiscreteFunction]==True,
(* Weights *)
If[
wflag==0,
Do[\[Lambda][i]=2HH Exp[x[i]^2]/(D[HermiteN[y],y])^2/.y->x[i],{i,1,n}];
weightsList=Table[\[Lambda][i],{i,1,n}];
Export[weightsName,weightsList];
Print[weightsName]];
wflag=1
];



];
(* End of Hermite Case *)


(* Beginning of Lengedre Case *)

If[type=="Legendre",
V[x_]:=potential/.variable->x;
t[y_]:=1/2 (a+b)+1/2 (b-a)y;
meshList=t[meshList];
Do[x[i]=meshList[[i]],{i,1,Length[meshList]}];

Do[For[i=1,i<j,i++,Subscript[T, i,j]=((-1)^(i+j) (a(x[i]+x[j])+b(x[i]+x[j])-2x[i]x[j]-2 a b))/(\[Sqrt]((x[i]-a)(b-x[i]))\[Sqrt]((x[j]-a)(b-x[j])) (x[i]-x[j])^2);Subscript[T, j,i]=Subscript[T, i,j]],{j,1,n}];

Do[Subscript[T, i,i]=((a-b)^2+n(n+1)(x[i]-a)(b-x[i]))/(3(a-x[i])^2 (b-x[i])^2),{i,1,n}];
Clear[i];


(* EigenFunctions *)
LegendreN[y_]:=(1/(b-a))^n (2n)!/(n!)^2 \!\(
\*UnderoverscriptBox[\(\[Product]\), \(k = 1\), \(n\)]\((y - x[k])\)\);
f[j_,y_]:=1/Sqrt[2] (-1)^(j+1) 1/(Sqrt[(x[j]-a)(b-x[j])]) LegendreN[y]/(y-x[j]);
\[Omega][y_]:=((y-a)(b-y))Sqrt[2/(b-a)];
If[OptionValue[DiscreteFunction]==True,
(* Weights *)
If[
wflag==0,
Do[\[Lambda][i]=2/((x[i]-a)(b-x[i])(D[LegendreN[y],y])^2)/.y->x[i],{i,1,n}];
weightsList=Table[\[Lambda][i],{i,1,n}];
Export[weightsName,weightsList];
Print[weightsName];
wflag=1
]
];

];






(* End of Legendre Case *)



(* Beginning of Laguerre Case *)
If[type=="Laguerre",
If[a<\[Infinity],
V[x_]:=potential/.variable->x+a;
V2[x_]:=EVpotential/.variable->x+a];
If[b<\[Infinity],
V[x_]:=potential/.variable->x+b;
V2[x_]:=EVpotential/.variable->x+b];
Do[x[i]=meshList[[i]],{i,1,Length[meshList]}];

Do[For[i=1,i<j,i++,Subscript[T, i,j]=(-1)^(i-j) (x[i]+x[j])/((x[i]x[j])^(1/2) (x[i]-x[j])^2);Subscript[T, j,i]=Subscript[T, i,j]],{j,1,n}];
Do[Subscript[T, i,i]=(4+(4n+2)x[i]-x[i]^2)/(12 x[i]^2),{i,1,n}];

Clear[i];

(* EigenFunctions *)
LaguerreN[y_]:= 1/n! \!\(
\*UnderoverscriptBox[\(\[Product]\), \(k = 1\), \(n\)]\((y - x[k])\)\);
f[j_,y_]:= (-1)^(j+n) x[j]^(-1/2) y LaguerreN[y]/(y-x[j]);
\[Omega][y_]:=Exp[-\[Sigma]((y)/2)];


If[OptionValue[DiscreteFunction]==True,
(* Weights *)
If[
wflag==0,
Do[\[Lambda][i]=Exp[x[i]]/(x[i](D[LaguerreN[y],y])^2)/.y->x[i],{i,1,n}];
weightsList=Table[\[Lambda][i],{i,1,n}];
Export[weightsName,weightsList];
Print[weightsName];
wflag=1]
];

];


(* End of Laguerre Case *)

min=0;
Quiet[If[Abs[
NMinimize[{V[x],a<=x<=b},x,WorkingPrecision->10][[1]]]<\[Infinity],
min=Rationalize[NMinimize[{V[x],a<=x<=b},x,WorkingPrecision->OptionValue[WorkingPrecision]][[1]],10^-OptionValue[WorkingPrecision]],
min=0]];

(* Matrix Potential Elements *)
Quiet[
Do[Subscript[V,i,i]=V[\[Sigma] h x[i]],{i,1,n}]];
U=DiagonalMatrix[Table[Subscript[V,i,i],{i,1,n}]];


(* Hamiltonian Matrix *)
H=Table[Table[1/(2 \[Mu] h^2) Subscript[T,l,m],{m,1,n}],{l,1,n}]+U;



(* Diagonalization *)

spectrum=Eigensystem[H-(min-PS)IdentityMatrix[n],-NLevels,Method->OptionValue[Method],Cubics->OptionValue[Cubics],Quartics->OptionValue[Quartics]];
energies=Table[spectrum[[1]][[NLevels-k]]+min-PS,{k,0,NLevels-1}];

Do[Do[\!\(\*OverscriptBox[
SubscriptBox[\(c\), \(j\)], \(k\)]\)=spectrum[[2]][[k]][[j]],{j,1,n}],{k,1,NLevels}];

(* Continuous Wave function *)

If[wflag==0,
If[OptionValue[CoefficientsOnly]==False,
wfunctions=Table[h^(-(1/2)) \[Omega][(variable+shift)/h]SetPrecision[Expand[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(
\*OverscriptBox[
SubscriptBox[\(c\), \(i\)], \(NLevels - k\)] f[i, 
\*FractionBox[\(\[Sigma] \((variable + shift)\)\), \(h\)]]\)\)],OptionValue[WorkingPrecision]],{k,0,NLevels-1}],
wfunctions=Table[Table[\!\(\*OverscriptBox[
SubscriptBox[\(c\), \(i\)], \(NLevels - k\)]\),{i,1,n}],{k,0,NLevels-1}]]
];

(* Discrete Wave function *)
If[wflag==1,
If[type=="Legendre",fac=(-1)^(n+1) Sqrt[2/(b-a)],fac=1];
wfunctions=Table[Table[{\[Sigma] h x[i]-h shift,fac/h^(1/2) \!\(\*OverscriptBox[
SubscriptBox[\(c\), \(i\)], \(NLevels - k\)]\)/Sqrt[\[Lambda][i]]},{i,1,n}],{k,0,NLevels-1}]];

If[SameQ[EVpotential,None]==True,
If[NLevels==1,Out={energies[[1]],wfunctions[[1]]},Out=Transpose[{energies,wfunctions}]]
];

(* Expectation Value *)
If[SameQ[EVpotential,None]==False,
EVlst=Table[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(
\*SuperscriptBox[\((
\*OverscriptBox[
SubscriptBox[\(c\), \(i\)], \(NLevels - k\)])\), \(2\)] V2[\[Sigma]\ h\ x[i]]\)\),{k,0,NLevels-1}];
If[NLevels==1,
Out={energies[[1]],wfunctions[[1]],EVlst[[1]]},
Out=Transpose[{energies,wfunctions,EVlst}]]
];
Out
]


(* ::Section::Closed:: *)
(*BuildMesh[]*)


(* Block for command BuildMesh *)

Options[BuildMesh]={WorkingPrecision->MachinePrecision,Weights->False};

BuildMesh[type_,n_,OptionsPattern[]]:=Module[{y,f,roots,dir,HH,HermiteN,LegendreN,LaguerreN,x,\[Lambda],weights},

(* Begin[Alerts] *)
BuildMesh::nnarg1 = "Value of option `1` is not a positive integer number larger than 1.";
BuildMesh::nnarg2 = "The mesh `1` does not coincide with Hermite, Legendre, or Laguerre.";
BuildMesh::nnarg3 = "Value of option `2` -> `1` is not a positive integer number larger than 0.";
BuildMesh::nnarg4 = "Value of option `2` -> `1` must be False or True.";

(* End[Alerts] *)


(* Begin[Abort] *)
If[type!="Laguerre"&&type!="Hermite"&&type!="Legendre",Message[BuildMesh::nnarg2, type];Abort[]];
If[SameQ[Simplify[Im[n]],0]==False\[Or]IntegerQ[n]==False\[Or]n<=1,Message[BuildMesh::nnarg1, n];Abort[]];
If[
NumericQ[OptionValue[WorkingPrecision]]==False, 
Message[BuildMesh::nnarg3, OptionValue[WorkingPrecision],WorkingPrecision];Abort[]];
If[
OptionValue[WorkingPrecision]!=MachinePrecision,
If[
IntegerQ[OptionValue[WorkingPrecision]]==False\[Or]OptionValue[WorkingPrecision]<1, 
Message[BuildMesh::nnarg3, OptionValue[WorkingPrecision],WorkingPrecision];Abort[]]];

If[MemberQ[{False,True},OptionValue[Weights]]==False,Message[BuildMesh::nnarg4, OptionValue[Weights],Weights];Abort[]]

(* End[Abort]*)


If[type=="Hermite",f=HermiteH[n,y]];
If[type=="Laguerre",f=LaguerreL[n,y]];
If[type=="Legendre",f=LegendreP[n,y]];


SetDirectory[NotebookDirectory[]];
dir=Directory[];

If[DirectoryQ["Meshes"]==False, CreateDirectory["Meshes"];
CreateDirectory["Meshes/Hermite/Weights/"];
CreateDirectory["Meshes/Laguerre/Weights/"];
CreateDirectory["Meshes/Legendre//Weights/"];
CreateDirectory["Meshes/Hermite/MeshPoints/"];
CreateDirectory["Meshes/Laguerre/MeshPoints/"];
CreateDirectory["Meshes/Legendre/MeshPoints/"]
];


roots=y/.NSolve[f==0, y, OptionValue[WorkingPrecision]];

Export[ToString[dir]<>"/Meshes/"<>ToString[type]<>"/MeshPoints/"<>ToString[type]<>"_"<>ToString[n]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>".dat",roots];
Print[ToString[type]<>"_"<>ToString[n]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>".dat"];



(* Quadrature's Weights *)
If[OptionValue[Weights]==True,
Do[x[i]=roots[[i]],{i,1,Length[roots]}];
(* Beginning  of Hermite Case *)
If[type=="Hermite",
(* Numerical Hermite Polynomials *)
HH=Sqrt[\[Pi]] 2^n n!;
HermiteN[y_]:=2^n \!\(
\*UnderoverscriptBox[\(\[Product]\), \(k = 1\), \(n\)]\((y - x[k])\)\);
Do[\[Lambda][i]=2HH Exp[x[i]^2]/(D[HermiteN[y],y])^2/.y->x[i],{i,1,n}]];
(* End of Hermite Case *)


(* Beginning of Lengedre Case *)
If[type=="Legendre",
(* Numerical Legendre Polynomials *)
LegendreN[y_]:=(1/2)^n (2n)!/(n!)^2 \!\(
\*UnderoverscriptBox[\(\[Product]\), \(k = 1\), \(n\)]\((y - x[k])\)\);
(* Quadrature's Weights *)

Do[\[Lambda][i]=2/((x[i]+1)(1-x[i])(D[LegendreN[y],y])^2)/.y->x[i],{i,1,n}]];
(* End of Legendre Case *)


(* Beginning of Laguerre Case *)
If[type=="Laguerre",
(* Numerical Hermite Laguerre *)
LaguerreN[y_]:= 1/n! \!\(
\*UnderoverscriptBox[\(\[Product]\), \(k = 1\), \(n\)]\((y - x[k])\)\);
(* Quadrature's Weights *)

Do[\[Lambda][i]=Exp[x[i]]/(x[i](D[LaguerreN[y],y])^2)/.y->x[i],{i,1,n}]];
(* End of Laguerre Case *)

weights=Table[\[Lambda][i],{i,1,n}];
Export[ToString[dir]<>"/Meshes/"<>ToString[type]<>"/Weights/"<>ToString[type]<>"_"<>ToString[n]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>"_Weights.dat",weights];
Print[ToString[type]<>"_"<>ToString[n]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>"_Weights.dat"];


];




]






(* ::Section:: *)
(*AvailableMeshQ[]*)


(* Block for AvailableMeshQ *)


DefaultPrecision=10^1000;
Options[AvailableMeshQ]={WorkingPrecision->DefaultPrecision,Dimension->DefaultPrecision,PrintMesh->False,PrintDomain->False};

AvailableMeshQ[type_,OptionsPattern[]]:=Module[{lst,lstOut,lstOutMP,dir,meshName},

SetDirectory[NotebookDirectory[]];
dir=Directory[];

(* Begin[Alerts] *)
AvailableMeshQ::nnarg1 = "The argument `1` does not coincide with Hermite, Legendre, or Laguerre.";
AvailableMeshQ::nnarg2 = "Value of option `1` -> "<>ToString[OptionValue[Dimension]]<>"   is not a positive integer number larger than 1.";
AvailableMeshQ::nnarg3 = "Value of option `1` -> "<>ToString[OptionValue[WorkingPrecision]]<>" is not a positive integer number larger than 0.";
AvailableMeshQ::nnarg4 = "The value of option `1` must be either False or True.";
AvailableMeshQ::nnarg5 = "The value of option `1` must be specified.";
AvailableMeshQ::nnarg6 = "The value of options `2` and `1` cannot be simultaneously True.";
(* End[Alerts] *)

(* Begin[Abort] *)
If[type!="Laguerre"&&type!="Hermite"&&type!="Legendre",Message[AvailableMeshQ::nnarg1, type];Abort[]];
dir=ToString[dir]<>"/Meshes/"<>type<>"/MeshPoints/";
SetDirectory[dir];

If[NumericQ[OptionValue[WorkingPrecision]]==False,Message[AvailableMeshQ::nnarg3, WorkingPrecision];Abort[]];
If[OptionValue[WorkingPrecision]!=MachinePrecision,
If[OptionValue[WorkingPrecision]<1\[Or]IntegerQ[OptionValue[WorkingPrecision]]==False,
Message[AvailableMeshQ::nnarg3, WorkingPrecision];Abort[]]
];
If[NumericQ[OptionValue[Dimension]]==False,Message[AvailableMeshQ::nnarg2, Dimension];Abort[]];
If[SameQ[Simplify[Im[OptionValue[Dimension]]],0]==False\[Or]OptionValue[Dimension]<=1\[Or]IntegerQ[OptionValue[Dimension]]==False,
Message[AvailableMeshQ::nnarg2, Dimension];Abort[]];
If[MemberQ[{False,True},OptionValue[PrintMesh]]==False,Message[AvailableMeshQ::nnarg4, PrintMesh];Abort[]];
If[MemberQ[{False,True},OptionValue[PrintDomain]]==False,Message[AvailableMeshQ::nnarg4, PrintDomain];Abort[]];
If[OptionValue[PrintMesh]==True&&OptionValue[PrintDomain]==True, Message[AvailableMeshQ::nnarg6, PrintMesh, PrintDomain];Abort[]];

(* End[Abort] *)

meshName=ToString[type]<>"_"<>ToString[OptionValue[Dimension]]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>".dat";

If[OptionValue[PrintMesh]==False,
If[FileExistsQ[meshName]&&OptionValue[PrintDomain]==False,Return[True,Module]]];

If[OptionValue[PrintDomain]==False,
If[FileExistsQ[meshName]&&OptionValue[PrintMesh]==False,Return[True,Module]]];

If[OptionValue[Dimension]!=DefaultPrecision&&OptionValue[WorkingPrecision]!=DefaultPrecision, 
If[FileExistsQ[ToString[type]<>"_"<>ToString[OptionValue[Dimension]]<>"_WP_"<>ToString[OptionValue[WorkingPrecision]]<>".dat"]==False,
Return[False]]
];

If[OptionValue[PrintMesh]==True,
If[OptionValue[WorkingPrecision]==DefaultPrecision,Message[AvailableMeshQ::nnarg5, WorkingPrecision];Abort[]];
If[OptionValue[Dimension]==DefaultPrecision,Message[AvailableMeshQ::nnarg5, Dimension];Abort[]]
];


If[OptionValue[PrintMesh]==True,
Return[Import[meshName,"List"],Module]
];


If[OptionValue[PrintDomain]==True,
If[OptionValue[WorkingPrecision]==DefaultPrecision,Message[AvailableMeshQ::nnarg5, WorkingPrecision];Abort[]];
If[OptionValue[Dimension]==DefaultPrecision,Message[AvailableMeshQ::nnarg5, Dimension];Abort[]]
];

If[OptionValue[PrintDomain]==True,
Return[{N[Import[meshName,"List"][[1]]],N[Import[meshName,"List"][[-1]]]},Module]
];


lstOut={};
lstOutMP={};

lst=FileNames["*.dat"];
ResetDirectory[];
dir=Directory[];
dir=ToString[dir]<>"/Meshes/"<>type<>"/Weights/";
SetDirectory[dir];
Do[
If[StringSplit[lst[[i]],{"_","."}][[4]]!="MachinePrecision",
AppendTo[lstOut,{ToExpression[StringSplit[lst[[i]],{"_","."}][[2]]],ToExpression[StringSplit[lst[[i]],{"_","."}][[4]]],
If[
FileExistsQ[
ToString[type]<>"_"<>ToString[ToExpression[StringSplit[lst[[i]],{"_","."}][[2]]]]<>"_WP_"<>ToString[ToExpression[StringSplit[lst[[i]],{"_","."}][[4]]]]<>"_Weights.dat"],
"Yes","No"]
}]
,AppendTo[lstOutMP,{ToExpression[StringSplit[lst[[i]],{"_","."}][[2]]],ToExpression[StringSplit[lst[[i]],{"_","."}][[4]]],
If[
FileExistsQ[
ToString[type]<>"_"<>ToString[ToExpression[StringSplit[lst[[i]],{"_","."}][[2]]]]<>"_WP_"<>ToString[ToExpression[StringSplit[lst[[i]],{"_","."}][[4]]]]<>"_Weights.dat"],
"Yes","No"]


}]
]
,{i,1,Length[lst]}];

lstOut=Flatten[Sort[ #, #1[[2]] < #2[[2]] &] & /@ SplitBy[Sort[lstOut, #1[[1]] < #2[[1]] &], First], 1];
lstOutMP=Sort[lstOutMP, #1[[1]] < #2[[1]] &];

lstOut=Join[lstOut,lstOutMP];  
If[OptionValue[WorkingPrecision]==DefaultPrecision&&OptionValue[Dimension]!=DefaultPrecision,
lstOut=Select[lstOut,First/*EqualTo[OptionValue[Dimension]]]
];

If[OptionValue[WorkingPrecision]!=DefaultPrecision&&OptionValue[Dimension]==DefaultPrecision,
lstOut=Select[lstOut, #[[2]]==OptionValue[WorkingPrecision]&]
];

If[Length[lstOut]==0,Return[False]];



TableForm[lstOut,TableHeadings->{None,{"Dimension","WorkingPrecision","Weights"}},TableAlignments -> Center]




]


(* ::Section:: *)
(*Protected Names*)


SetAttributes[{Scaling,Dimension,Mass,PotentialShift,ExpectationValue,Weights,DiscreteFunction,PrintMesh,PrintDomain},Protected]


End[];


 EndPackage[];


(* ::Subsection:: *)
(**)
