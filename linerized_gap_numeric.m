(* ::Package:: *)

ClearAll["`*"];


LaunchKernels[6];


IR[\[Eta]_]=Floor[-Log[Exp[\[Eta]]-1]/\[Eta]]+1;offset[\[Eta]_]=IR[\[Eta]]-Floor[Log[IR[\[Eta]]]/\[Eta]]-1;
e[n_,\[Eta]_]:=If[n<IR[\[Eta]],n,Floor[Exp[\[Eta]*(n-offset[\[Eta]])]]];
g=1;\[Xi][T_]=2/3 (g^2/(2\[Pi] T))^(1/3);
\[Xi]1[T_]=(3 \[Pi]^(1/3) (T/g^2)^(1/3))/2^(2/3);
\[CapitalSigma]T[T_]=(3 g^(1/3))/((2 \[Pi])^(2/3) T^(1/6));

K1=Compile[{{m,_Real},{n,_Real},{T,_Real},{\[Eta],_Real}},If[m==n,(e[n+1,\[Eta]]-e[n,\[Eta]])/(\[Xi]1[T](2e[n,\[Eta]]-1)+1000*\[CapitalSigma]T[T]+2(Zeta[1/3]-Zeta[1/3,e[n,\[Eta]]]) ) 1/(2e[n,\[Eta]]-1)^(1/3),Sqrt[e[n+1,\[Eta]]-e[n,\[Eta]]]/Abs[\[Xi]1[T](2e[n,\[Eta]]-1)+1000*\[CapitalSigma]T[T]+2(Zeta[1/3]-Zeta[1/3,e[n,\[Eta]]])]^(1/2) (1/Abs[e[m,\[Eta]]-e[n,\[Eta]]]^(1/3)+1/(e[m,\[Eta]]+e[n,\[Eta]]-1)^(1/3))Sqrt[e[m+1,\[Eta]]-e[m,\[Eta]]]/Abs[\[Xi]1[T](2e[m,\[Eta]]-1)+1000*\[CapitalSigma]T[T]+2(Zeta[1/3]-Zeta[1/3,e[m,\[Eta]]])]^(1/2) ]];
maxEig[T_,Nmax_,\[Eta]_]:=Block[{vin,Kf,Ktemp,KT},
vin=Table[{i,j},{i,1,Nmax},{j,1,i}];
vin=Flatten[vin,1];
Kf[x_]:=Kf[x]=K1[x[[1]],x[[2]],T,\[Eta]];
Ktemp=ParallelMap[Kf,vin];
KT=Table[0,{i,1,Nmax},{j,1,Nmax}];
Do[KT[[i,j]]=KT[[j,i]]=Ktemp[[1/2 (i-1)i+j]],{i,1.,Nmax},{j,1.,i}];
Eigensystem[KT,1,Method->{Arnoldi,MaxIterations->10^2}]];


(* T=1*)

tT1Nmax2000eta005=maxEig[1. ,2000,0.05];
Export["tT1Nmax2000eta005.dat",tT1Nmax2000eta005];
Clear[tT1Nmax2000eta005];

tT1Nmax2500eta005=maxEig[1. ,2500,0.05];
Export["tT1Nmax2500eta005.dat",tT1Nmax2500eta005];
Clear[tT1Nmax2500eta005];

tT1Nmax1500eta007=maxEig[1. ,1500,0.07];
Export["tT1Nmax1500eta007.dat",tT1Nmax1500eta007];
Clear[tT1Nmax1500eta007];

tT1Nmax2000eta007=maxEig[1. ,2000,0.07];
Export["tT1Nmax2000eta007.dat",tT1Nmax2000eta007];
Clear[tT1Nmax2000eta007];

tT1Nmax1500eta009=maxEig[1. ,1500,0.09];
Export["tT1Nmax1500eta009.dat",tT1Nmax1500eta009];
Clear[tT1Nmax1500eta009];

tT1Nmax2000eta009=maxEig[1. ,2000,0.09];
Export["tT1Nmax2000eta009.dat",tT1Nmax2000eta009];
Clear[tT1Nmax2000eta009];

tT1Nmax1500eta01=maxEig[1. ,1500,0.1];
Export["tT1Nmax1500eta01.dat",tT1Nmax1500eta01];
Clear[tT1Nmax1500eta01];

tT1Nmax2000eta01=maxEig[1. ,2000,0.1];
Export["tT1Nmax2000eta01.dat",tT1Nmax2000eta01];
Clear[tT1Nmax2000eta01];


(* T=10^-10*)

tTE10Nmax2000eta005=maxEig[1. 10^-10,2000,0.05];
Export["tTE10Nmax2000eta005.dat",tTE10Nmax2000eta005];
Clear[tTE10Nmax2000eta005];

tTE10Nmax2500eta005=maxEig[1. 10^-10,2500,0.05];
Export["tTE10Nmax2500eta005.dat",tTE10Nmax2500eta005];
Clear[tTE10Nmax2500eta005];

tTE10Nmax2000eta007=maxEig[1. 10^-10,2000,0.07];
Export["tTE10Nmax2000eta007.dat",tTE10Nmax2000eta007];
Clear[ttTE10Nmax2000eta007];

tTE10Nmax2500eta007=maxEig[1. 10^-10,2500,0.07];
Export["tTE10Nmax2500eta007.dat",tTE10Nmax2500eta007];
Clear[ttTE10Nmax2500eta007];

tTE10Nmax1500eta009=maxEig[1. 10^-10,1500,0.09];
Export["tTE10Nmax1500eta009.dat",tTE10Nmax1500eta009];
Clear[tTE10Nmax1500eta009];

tTE10Nmax2000eta009=maxEig[1. 10^-10,2000,0.09];
Export["tTE10Nmax2000eta009.dat",tTE10Nmax2000eta009];
Clear[tTE10Nmax2000eta009];

tTE10Nmax1500eta01=maxEig[1. 10^-10,1500,0.1];
Export["tTE10Nmax1500eta01.dat",tTE10Nmax1500eta01];
Clear[tTE10Nmax1500eta01];

tTE10Nmax2000eta01=maxEig[1. 10^-10,2000,0.1];
Export["tTE10Nmax2000eta01.dat",tTE10Nmax2000eta01];
Clear[ttTE10Nmax2000eta01];


(* T=10^-20*)

tTE20Nmax2000eta005=maxEig[1. 10^-20,2000,0.05];
Export["tTE20Nmax2000eta005.dat",tTE20Nmax2000eta005];
Clear[tTE20Nmax2000eta005];

tTE20Nmax2500eta005=maxEig[1. 10^-20,2500,0.05];
Export["tTE20Nmax2500eta005.dat",tTE20Nmax2500eta005];
Clear[tTE20Nmax2500eta005];

tTE20Nmax2000eta007=maxEig[1. 10^-20,2000,0.07];
Export["tTE20Nmax2000eta007.dat",tTE20Nmax2000eta007];
Clear[tTE20Nmax2000eta007];

tTE20Nmax2500eta007=maxEig[1. 10^-20,2500,0.07];
Export["tTE20Nmax2500eta007.dat",tTE20Nmax2500eta007];
Clear[tTE20Nmax2500eta007];

tTE20Nmax1500eta009=maxEig[1. 10^-20,1500,0.09];
Export["tTE20Nmax1500eta009.dat",tTE20Nmax1500eta009];
Clear[tTE20Nmax1500eta009];

tTE20Nmax2000eta009=maxEig[1. 10^-20,2000,0.09];
Export["tTE20Nmax2000eta009.dat",tTE20Nmax2000eta009];
Clear[tTE20Nmax2000eta009];

tTE20Nmax1500eta01=maxEig[1. 10^-20,1500,0.1];
Export["tTE20Nmax1500eta01.dat",tTE20Nmax1500eta01];
Clear[tTE20Nmax1500eta01];

tTE20Nmax2000eta01=maxEig[1. 10^-20,2000,0.1];
Export["tTE20Nmax2000eta01.dat",tTE20Nmax2000eta01];
Clear[tTE20Nmax2000eta01];


(* T=10^-30*)

tTE30Nmax2000eta005=maxEig[1. 10^-30,2000,0.05];
Export["tTE30Nmax2000eta005.dat",tTE30Nmax2000eta005];
Clear[tTE30Nmax2000eta005];

tTE30Nmax2500eta005=maxEig[1. 10^-30,2500,0.05];
Export["tTE30Nmax2500eta005.dat",tTE30Nmax2500eta005];
Clear[tTE30Nmax2500eta005];

tTE30Nmax2000eta007=maxEig[1. 10^-30,2000,0.07];
Export["tTE30Nmax2000eta007.dat",tTE30Nmax2000eta007];
Clear[tTE30Nmax2000eta007];

tTE30Nmax2500eta007=maxEig[1. 10^-30,2500,0.07];
Export["tTE30Nmax2500eta007.dat",tTE30Nmax2500eta007];
Clear[tTE30Nmax2500eta007];

tTE30Nmax2000eta009=maxEig[1. 10^-30,2000,0.09];
Export["tTE30Nmax2000eta009.dat",tTE30Nmax2000eta009];
Clear[tTE30Nmax2000eta009];

tTE30Nmax2500eta009=maxEig[1. 10^-30,2500,0.09];
Export["tTE30Nmax2500eta009.dat",tTE30Nmax2500eta009];
Clear[tTE30Nmax2500eta009];

tTE30Nmax1500eta01=maxEig[1. 10^-30,1500,0.1];
Export["tTE30Nmax1500eta01.dat",tTE30Nmax1500eta01];
Clear[tTE30Nmax1500eta01];

tTE30Nmax2000eta01=maxEig[1. 10^-30,2000,0.1];
Export["tTE30Nmax2000eta01.dat",tTE30Nmax2000eta01];
Clear[tTE30Nmax2000eta01];


(* T=10^-40*)

tTE40Nmax2000eta005=maxEig[1. 10^-40,2000,0.05];
Export["tTE40Nmax2000eta005.dat",tTE40Nmax2000eta005];
Clear[tTE40Nmax2000eta005];

tTE40Nmax2500eta005=maxEig[1. 10^-40,2500,0.05];
Export["tTE40Nmax2500eta005.dat",tTE40Nmax2500eta005];
Clear[tTE40Nmax2500eta005];

tTE40Nmax2000eta007=maxEig[1. 10^-40,2000,0.07];
Export["tTE40Nmax2000eta007.dat",tTE40Nmax2000eta007];
Clear[tTE40Nmax2000eta007];

tTE40Nmax2500eta007=maxEig[1. 10^-40,2500,0.07];
Export["tTE40Nmax2500eta007.dat",tTE40Nmax2500eta007];
Clear[tTE40Nmax2500eta007];

tTE40Nmax2000eta009=maxEig[1. 10^-40,2000,0.09];
Export["tTE40Nmax2000eta009.dat",tTE40Nmax2000eta009];
Clear[tTE40Nmax2000eta009];

tTE40Nmax2500eta009=maxEig[1. 10^-40,2500,0.09];
Export["tTE40Nmax2500eta009.dat",tTE40Nmax2500eta009];
Clear[tTE40Nmax2500eta009];

tTE40Nmax1500eta01=maxEig[1. 10^-40,1500,0.1];
Export["tTE40Nmax1500eta01.dat",tTE40Nmax1500eta01];
Clear[tTE40Nmax1500eta01];

tTE40Nmax2000eta01=maxEig[1. 10^-40,2000,0.1];
Export["tTE40Nmax2000eta01.dat",tTE40Nmax2000eta01];
Clear[tTE40Nmax2000eta01];


(* T=10^-50*)

tTE50Nmax2500eta005=maxEig[1. 10^-50,2500,0.05]
Export["tTE50Nmax2500eta005.dat",tTE50Nmax2500eta005];
Clear[tTE50Nmax2500eta005];

tTE50Nmax3000eta005=maxEig[1. 10^-50,3000,0.05];
Export["tTE50Nmax3000eta005.dat",tTE50Nmax3000eta005];
Clear[tTE50Nmax3000eta005];

tTE50Nmax2500eta007=maxEig[1. 10^-50,2500,0.07];
Export["tTE50Nmax2500eta007.dat",tTE50Nmax2500eta007];
Clear[tTE50Nmax2500eta007];

tTE50Nmax3000eta007=maxEig[1. 10^-50,3000,0.07];
Export["tTE50Nmax3000eta007.dat",tTE50Nmax3000eta007];
Clear[tTE50Nmax3000eta007];

tTE50Nmax2000eta009=maxEig[1. 10^-50,2000,0.09];
Export["tTE50Nmax2000eta009.dat",tTE50Nmax2000eta009];
Clear[tTE50Nmax2000eta009];

tTE50Nmax2500eta009=maxEig[1. 10^-50,2500,0.09];
Export["tTE50Nmax2500eta009.dat",tTE50Nmax2500eta009];
Clear[tTE50Nmax2500eta009];

tTE50Nmax1500eta01=maxEig[1. 10^-50,1500,0.1];
Export["tTE50Nmax1500eta01.dat",tTE50Nmax1500eta01];
Clear[tTE50Nmax1500eta01];

tTE50Nmax2000eta01=maxEig[1. 10^-50,2000,0.1];
Export["tTE50Nmax2000eta01.dat",tTE50Nmax2000eta01];
Clear[tTE50Nmax2000eta01];


(* T=10^-60*)

tTE60Nmax2500eta005=maxEig[1. 10^-60,2500,0.05];
Export["tTE60Nmax2500eta005.dat",tTE60Nmax2500eta005];
Clear[tTE60Nmax2500eta005];

tTE60Nmax2500eta005=maxEig[1. 10^-60,3000,0.05];
Export["tTE60Nmax2500eta005.dat",tTE60Nmax2500eta005];
Clear[tTE60Nmax2500eta005];

tTE60Nmax3500eta005=maxEig[1. 10^-60,3500,0.05];
Export["tTE60Nmax3500eta005.dat",tTE60Nmax3500eta005];
Clear[tTE60Nmax3500eta005];

tTE60Nmax2500eta007=maxEig[1. 10^-60,2500,0.07];
Export["tTE60Nmax2500eta007.dat",tTE60Nmax2500eta007];
Clear[tTE60Nmax2500eta007];

tTE60Nmax3000eta007=maxEig[1. 10^-60,3000,0.07];
Export["tTE60Nmax3000eta007.dat",tTE60Nmax3000eta007];
Clear[tTE60Nmax3000eta007];

tTE60Nmax2000eta009=maxEig[1. 10^-60,2000,0.09];
Export["tTE60Nmax2000eta009.dat",tTE60Nmax2000eta009];
Clear[tTE60Nmax2000eta009];

tTE60Nmax2500eta009=maxEig[1. 10^-60,2500,0.09];
Export["tTE60Nmax2500eta009.dat",tTE60Nmax2500eta009];
Clear[tTE60Nmax2500eta009];

tTE60Nmax1500eta01=maxEig[1. 10^-60,1500,0.1];
Export["tTE60Nmax1500eta01.dat",tTE60Nmax1500eta01];
Clear[tTE60Nmax1500eta01];

tTE60Nmax2000eta01=maxEig[1. 10^-60,2000,0.1];
Export["tTE60Nmax2000eta01.dat",tTE60Nmax2000eta01];
Clear[tTE60Nmax2000eta01];

tTE60Nmax1000eta02=maxEig[1. 10^-60,1000,0.2];
Export["tTE60Nmax1000eta02.dat",tTE60Nmax1000eta02];
Clear[tTE60Nmax1000eta02];

tTE60Nmax1500eta02=maxEig[1. 10^-60,1500,0.2];
Export["tTE60Nmax1500eta02.dat",tTE60Nmax1500eta02];
Clear[tTE60Nmax1500eta02];


(* T=10^-70*)

tTE70Nmax3000eta005=maxEig[1. 10^-70,3000,0.05];
Export["tTE70Nmax3000eta005.dat",tTE70Nmax3000eta005];
Clear[tTE70Nmax3000eta005];

tTE70Nmax3500eta005=maxEig[1. 10^-70,3500,0.05];
Export["tTE70Nmax3500eta005.dat",tTE70Nmax3500eta005];
Clear[tTE70Nmax3500eta005];

tTE70Nmax2500eta007=maxEig[1. 10^-70,2500,0.07];
Export["tTE70Nmax2500eta007.dat",tTE70Nmax2500eta007];
Clear[tTE70Nmax2500eta007];

tTE70Nmax3000eta007=maxEig[1. 10^-70,3000,0.07];
Export["tTE70Nmax3000eta007.dat",tTE70Nmax3000eta007];
Clear[tTE70Nmax3000eta007];

tTE70Nmax2000eta009=maxEig[1. 10^-70,2000,0.09];
Export["tTE70Nmax2000eta009.dat",tTE70Nmax2000eta009];
Clear[tTE70Nmax2000eta009];

tTE70Nmax2500eta009=maxEig[1. 10^-70,2500,0.09];
Export["tTE70Nmax2500eta009.dat",tTE70Nmax2500eta009];
Clear[tTE70Nmax2500eta009];

tTE70Nmax1500eta01=maxEig[1. 10^-70,1500,0.1];
Export["tTE70Nmax1500eta01.dat",tTE70Nmax1500eta01];
Clear[tTE70Nmax1500eta01];

tTE70Nmax2000eta01=maxEig[1. 10^-70,2000,0.1];
Export["tTE70Nmax2000eta01.dat",tTE70Nmax2000eta01];
Clear[tTE70Nmax2000eta01];


(* T=10^-80*)

tTE80Nmax2500eta007=maxEig[1. 10^-80,2500,0.07];
Export["tTE80Nmax2500eta007.dat",tTE80Nmax2500eta007];
Clear[tTE80Nmax2500eta007];

tTE80Nmax3000eta007=maxEig[1. 10^-80,3000,0.07];
Export["tTE80Nmax3000eta007.dat",tTE80Nmax3000eta007];
Clear[tTE80Nmax3000eta007];

tTE80Nmax2000eta009=maxEig[1. 10^-80,2000,0.09];
Export["tTE80Nmax2000eta009.dat",tTE80Nmax2000eta009];
Clear[tTE80Nmax2000eta009];

tTE80Nmax2500eta009=maxEig[1. 10^-80,2500,0.09];
Export["tTE80Nmax2500eta009.dat",tTE80Nmax2500eta009];
Clear[tTE80Nmax2500eta009];

tTE80Nmax2000eta01=maxEig[1. 10^-80,2000,0.1];
Export["tTE80Nmax2000eta01.dat",tTE80Nmax2000eta01];
Clear[tTE80Nmax2000eta01];

tTE80Nmax2500eta01=maxEig[1. 10^-80,2500,0.1];
Export["tTE80Nmax2500eta01.dat",tTE80Nmax2500eta01];
Clear[tTE80Nmax2500eta01];


(* T=10^-90*)

tTE90Nmax2500eta007=maxEig[1. 10^-90,2500,0.07];
Export["tTE90Nmax2500eta007.dat",tTE90Nmax2500eta007];
Clear[tTE90Nmax2500eta007];

tTE90Nmax3000eta007=maxEig[1. 10^-90,3000,0.07];
Export["tTE90Nmax3000eta007.dat",tTE90Nmax3000eta007];
Clear[tTE90Nmax3000eta007];

tTE90Nmax2000eta009=maxEig[1. 10^-90,2000,0.09];
Export["tTE90Nmax2000eta009.dat",tTE90Nmax2000eta009];
Clear[tTE90Nmax2000eta009];

tTE90Nmax2500eta009=maxEig[1. 10^-90,2500,0.09];
Export["tTE90Nmax2500eta009.dat",tTE90Nmax2500eta009];
Clear[tTE90Nmax2500eta009];

tTE90Nmax2000eta01=maxEig[1. 10^-90,2000,0.1];
Export["tTE90Nmax2000eta01.dat",tTE90Nmax2000eta01];
Clear[tTE90Nmax2000eta01];

tTE90Nmax2500eta01=maxEig[1. 10^-90,2500,0.1];
Export["tTE90Nmax2500eta01.dat",tTE90Nmax2500eta01];
Clear[tTE90Nmax2500eta01];

tTE90Nmax2000eta011=maxEig[1. 10^-90,2000,0.11];
Export["tTE90Nmax2000eta011.dat",tTE90Nmax2000eta011];
Clear[tTE90Nmax2000eta011];

tTE90Nmax2500eta011=maxEig[1. 10^-90,2500,0.11];
Export["tTE90Nmax2500eta011.dat",tTE90Nmax2500eta011];
Clear[tTE90Nmax2500eta011];


(* T=10^-100*)

tTE100Nmax3000eta007=maxEig[1. 10^-100,3000,0.07];
Export["tTE100Nmax3000eta007.dat",tTE100Nmax3000eta007];
Clear[tTE100Nmax3000eta007];

tTE100Nmax3500eta007=maxEig[1. 10^-100,3500,0.07];
Export["tTE100Nmax3500eta007.dat",tTE100Nmax3500eta007];
Clear[tTE100Nmax3500eta007];

tTE100Nmax2500eta009=maxEig[1. 10^-100,2500,0.09];
Export["tTE100Nmax2500eta009.dat",tTE100Nmax2500eta009];
Clear[tTE100Nmax2500eta009];

tTE100Nmax3000eta009=maxEig[1. 10^-100,3000,0.09];
Export["tTE100Nmax3000eta009.dat",tTE100Nmax3000eta009];
Clear[tTE100Nmax3000eta009];

tTE100Nmax2000eta01=maxEig[1. 10^-100,2000,0.1];
Export["tTE100Nmax2000eta01.dat",tTE100Nmax2000eta01];
Clear[tTE100Nmax2000eta01];

tTE100Nmax2500eta01=maxEig[1. 10^-100,2500,0.1];
Export["tTE100Nmax2500eta01.dat",tTE100Nmax2500eta01];
Clear[tTE100Nmax2500eta01];

tTE100Nmax2000eta011=maxEig[1. 10^-100,2000,0.11];
Export["tTE100Nmax2000eta011.dat",tTE100Nmax2000eta011];
Clear[tTE100Nmax2000eta011];

tTE100Nmax2500eta011=maxEig[1. 10^-100,2500,0.11];
Export["tTE100Nmax2500eta011.dat",tTE100Nmax2500eta011];
Clear[tTE100Nmax2500eta011];


(* T=10^-120*)

tTE120Nmax3000eta009=maxEig[1. 10^-120,3000,0.09];
Export["tTE120Nmax3000eta009.dat",tTE120Nmax3000eta009];
Clear[tTE120Nmax3000eta009];

tTE120Nmax3500eta009=maxEig[1. 10^-120,3500,0.09];
Export["tTE120Nmax3500eta009.dat",tTE120Nmax3500eta009];
Clear[tTE120Nmax3500eta009];

tTE120Nmax3000eta011=maxEig[1. 10^-120,3000,0.11];
Export["tTE120Nmax3000eta011.dat",tTE120Nmax3000eta011];
Clear[tTE120Nmax3000eta011];

tTE120Nmax3500eta011=maxEig[1. 10^-120,3500,0.11];
Export["tTE120Nmax3500eta011.dat",tTE120Nmax3500eta011];
Clear[tTE120Nmax3500eta011];

tTE120Nmax2500eta013=maxEig[1. 10^-120,2500,0.13];
Export["tTE120Nmax2500eta013.dat",tTE120Nmax2500eta013];
Clear[tTE120Nmax2500eta013];

tTE120Nmax3000eta013=maxEig[1. 10^-120,3000,0.13];
Export["tTE120Nmax3000eta013.dat",tTE120Nmax3000eta013];
Clear[tTE120Nmax3000eta013];



(* T=10^-140*)

tTE140Nmax3000eta009=maxEig[1. 10^-140,3000,0.09];
Export["tTE140Nmax3000eta009.dat",tTE140Nmax3000eta009];
Clear[tTE140Nmax3000eta009];

tTE140Nmax3500eta009=maxEig[1. 10^-140,3500,0.09];
Export["tTE140Nmax3500eta009.dat",tTE140Nmax3500eta009];
Clear[tTE140Nmax3500eta009];

tTE140Nmax3000eta011=maxEig[1. 10^-140,3000,0.11];
Export["tTE140Nmax3000eta011.dat",tTE140Nmax3000eta011];
Clear[ttTE140Nmax3000eta011];

tTE140Nmax3500eta011=maxEig[1. 10^-140,3500,0.11];
Export["tTE140Nmax3500eta011.dat",tTE140Nmax3500eta011];
Clear[tTE140Nmax3500eta011];

tTE140Nmax2500eta013=maxEig[1. 10^-140,2500,0.13];
Export["tTE140Nmax2500eta013.dat",tTE140Nmax2500eta013];
Clear[tTE140Nmax2500eta013];

tTE140Nmax3000eta013=maxEig[1. 10^-140,3000,0.13];
Export["tTE140Nmax3000eta013.dat",tTE140Nmax3000eta013];
Clear[tTE140Nmax3000eta013];


(* T=10^-160*)

tTE160Nmax3000eta009=maxEig[1. 10^-160,3000,0.09];
Export["tTE160Nmax3000eta009.dat",tTE160Nmax3000eta009];
Clear[tTE160Nmax3000eta009];

tTE160Nmax3500eta009=maxEig[1. 10^-160,3500,0.09];
Export["tTE160Nmax3500eta009.dat",tTE160Nmax3500eta009];
Clear[tTE160Nmax3500eta009];

tTE160Nmax3000eta011=maxEig[1. 10^-160,3000,0.11];
Export["tTE160Nmax3000eta011.dat",tTE160Nmax3000eta011];
Clear[tTE160Nmax3000eta011];

tTE160Nmax3500eta011=maxEig[1. 10^-160,3500,0.11];
Export["tTE160Nmax3500eta011.dat",tTE160Nmax3500eta011];
Clear[tTE160Nmax3500eta011];

tTE160Nmax2500eta013=maxEig[1. 10^-160,2500,0.13];
Export["tTE160Nmax2500eta013.dat",tTE160Nmax2500eta013];
Clear[tTE160Nmax2500eta013];

tTE160Nmax3000eta013=maxEig[1. 10^-160,3000,0.13];
Export["tTE160Nmax3000eta013.dat",tTE160Nmax3000eta013];
Clear[tTE160Nmax3000eta013];


(* T=10^-180*)

tTE180Nmax3000eta009=maxEig[1. 10^-180,3000,0.09];
Export["tTE180Nmax3000eta009.dat",tTE180Nmax3000eta009];
Clear[tTE180Nmax3000eta009];

tTE180Nmax3500eta009=maxEig[1. 10^-180,3500,0.09];
Export["tTE180Nmax3500eta009.dat",tTE180Nmax3500eta009];
Clear[tTE180Nmax3500eta009];

tTE180Nmax3000eta011=maxEig[1. 10^-180,3000,0.11];
Export["tTE180Nmax3000eta011.dat",tTE180Nmax3000eta011];
Clear[tTE180Nmax3000eta011];

tTE180Nmax3500eta011=maxEig[1. 10^-180,3500,0.11];
Export["tTE180Nmax3500eta011.dat",tTE180Nmax3500eta011];
Clear[tTE180Nmax3500eta011];

tTE180Nmax2500eta013=maxEig[1. 10^-180,2500,0.13];
Export["tTE180Nmax2500eta013.dat",tTE180Nmax2500eta013];
Clear[tTE180Nmax2500eta013];

tTE180Nmax3000eta013=maxEig[1. 10^-180,3000,0.13];
Export["tTE180Nmax3000eta013.dat",tTE180Nmax3000eta013];
Clear[tTE180Nmax3000eta013];


(* T=10^-200*)

tTE200Nmax3000eta009=maxEig[1. 10^-200,3000,0.09];
Export["tTE200Nmax3000eta009.dat",tTE200Nmax3000eta009];
Clear[tTE200Nmax3000eta009];

tTE200Nmax3500eta009=maxEig[1. 10^-200,3500,0.09];
Export["tTE200Nmax3500eta009.dat",tTE200Nmax3500eta009];
Clear[tTE200Nmax3500eta009];

tTE200Nmax3000eta011=maxEig[1. 10^-200,3000,0.11];
Export["tTE200Nmax3000eta011.dat",tTE200Nmax3000eta011];
Clear[tTE200Nmax3000eta011];

tTE200Nmax3500eta011=maxEig[1. 10^-200,3500,0.11];
Export["tTE200Nmax3500eta011.dat",tTE200Nmax3500eta011];
Clear[tTE200Nmax3500eta011];

tTE200Nmax2500eta013=maxEig[1. 10^-200,2500,0.13];
Export["tTE200Nmax2500eta013.dat",tTE200Nmax2500eta013];
Clear[tTE200Nmax2500eta013];

tTE200Nmax3000eta013=maxEig[1. 10^-200,3000,0.13];
Export["tTE200Nmax3000eta013.dat",tTE200Nmax3000eta013];
Clear[tTE200Nmax3000eta013];
