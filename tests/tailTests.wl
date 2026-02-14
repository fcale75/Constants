root=DirectoryName[$InputFileName];
Get[FileNameJoin[{root,"..","wl","Constants.wl"}]];
Get[FileNameJoin[{root,"..","wl","Newton.wl"}]];
Get[FileNameJoin[{root,"..","wl","TailDerivatives.wl"}]];
ConstantsSetPrecision[80];

P=4; Nmax=6; K=4;
a=ConstantArray[1/P,P];
eps=SetPrecision[10^-8,40];
obj[a_]:=Objective[a,Nmax,K];
numGrad=Table[(obj[a+eps UnitVector[P,j]]-obj[a-eps UnitVector[P,j]])/(2 eps),{j,1,P}];

rep=Constants`Newton`NewtonDerivatives[a,P,Nmax];
g=rep["gpart"]+Constants`TailDerivatives`TailGradient[a,Nmax,K];
h=rep["hess"]+Constants`TailDerivatives`TailHessian[a,Nmax,K];

tests={
VerificationTest[VectorQ[g,NumericQ]&&Length[g]==P],
VerificationTest[Norm[g-numGrad]<10^-25],
VerificationTest[Max[Abs[h-Transpose[h]]]<10^-40]
};
Print[TestReport[tests]];
