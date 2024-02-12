(* ::Package:: *)

ClearAll["`*"];


LaunchKernels[8];


tsol[x_]:=E^x Hypergeometric2F1[1/2-1/2 I Sqrt[-1+4 g1],1/2+1/2 I Sqrt[-1+4 g1],2,-E^x]+c MeijerG[{{},{1/2 (3-I Sqrt[-1+4 g1]),1/2 (3+I Sqrt[-1+4 g1])}},{{0,1},{}},-E^x]


csol[\[CapitalLambda]_,g1_]:=(E^-\[CapitalLambda] g1 Hypergeometric2F1[3/2-1/2 I Sqrt[-1+4 g1],3/2+1/2 I Sqrt[-1+4 g1],3,-E^-\[CapitalLambda]])/(2 (MeijerG[{{},{1/2 (1-I Sqrt[-1+4 g1]),1/2 (1+I Sqrt[-1+4 g1])}},{{0,0},{}},-E^-\[CapitalLambda]]-E^\[CapitalLambda] MeijerG[{{},{1/2 (3-I Sqrt[-1+4 g1]),1/2 (3+I Sqrt[-1+4 g1])}},{{0,1},{}},-E^-\[CapitalLambda]]));
sol[x_,g1_]:=E^x Hypergeometric2F1[1/2-1/2 I Sqrt[-1+4 g1],1/2+1/2 I Sqrt[-1+4 g1],2,-E^x]+c MeijerG[{{},{1/2 (3-I Sqrt[-1+4 g1]),1/2 (3+I Sqrt[-1+4 g1])}},{{0,1},{}},-E^x]/.c->csol[\[CapitalLambda],g1];


\[CapitalLambda]=10;
gs[m_,n_]:=1/4 1/((70+m)/10) (8+n/100);


nsol=ParallelTable[\[CapitalDelta]/.NDSolve[{\[CapitalDelta]'[x]-\[CapitalDelta]''[x]==gs[m,n]/(1+E^-x) \[CapitalDelta][x],\[CapitalDelta][-\[CapitalLambda]]==Re[sol[-\[CapitalLambda],gs[m,n]]],\[CapitalDelta]'[-\[CapitalLambda]]==Re[sol[-\[CapitalLambda],gs[m,n]]]},\[CapitalDelta],{x,-\[CapitalLambda],1000},WorkingPrecision->20][[1]],{n,1,40},{m,1,10}];


x0t[m_,n_]:=(2\[Pi])/Sqrt[4gs[m,n]-1]//N;


xt=ParallelTable[x/.FindRoot[nsol[[m]][[n]]'[x]==0,{x,x0t[m,n]},WorkingPrecision->20],{n,1,40},{m,1,10}];


\[Delta]=10^-20;
u[x_,xp_]:=1/Abs[E^(-3x)-E^(-3xp)]^(1/3)+1/Abs[E^(-3x)+E^(-3xp)]^(1/3);


(* w/ Sigma_T *)
Int[x_?NumericQ,m_?NumericQ,n_?NumericQ]:=1/((70+m)/10) NIntegrate[u[x,xp] nsol[[m]][[n]][xp]/(1+E^xp+E^(3xp-2xt[[m]][[n]])),{xp,-\[CapitalLambda],x-\[Delta]},WorkingPrecision->15,MaxRecursion->20]+1/((70+m)/10) NIntegrate[u[x,xp] nsol[[m]][[n]][xp]/(1+E^xp+E^(3xp-2xt[[m]][[n]])),{xp,x+\[Delta],xt[[m]][[n]]},WorkingPrecision->15,MaxRecursion->20];


(* w/out Sigma_T*)
Int2[x_?NumericQ,m_?NumericQ,n_?NumericQ]:=1/((70+m)/10) NIntegrate[u[x,xp] nsol[[m]][[n]][xp]/(1+E^xp),{xp,-\[CapitalLambda],x-\[Delta]},WorkingPrecision->15,MaxRecursion->20]+1/((70+m)/10) NIntegrate[u[x,xp] nsol[[m]][[n]][xp]/(1+E^xp),{xp,x+\[Delta],xt[[m]][[n]]},WorkingPrecision->15,MaxRecursion->20];


xi=-\[CapitalLambda];
Steps=40;
dx=ParallelTable[(xt[[m]][[n]]-xi)/(Steps-1),{n,1,40},{m,1,10}];
xlist=ParallelTable[Table[xi+dx[[m]][[n]](i-1),{i,1,Steps}],{n,1,40},{m,1,10}];


\[CapitalDelta]tab2=ParallelTable[Table[{xlist[[m]][[n]][[i]],Int2[xlist[[m]][[n]][[i]],m,n]},{i,1,Steps}],{n,15,40},{m,1,10}];
Export["\[CapitalDelta]tab2.dat",\[CapitalDelta]tab2];
Clear[\[CapitalDelta]tab2];


\[CapitalDelta]tab=ParallelTable[Table[{xlist[[m]][[n]][[i]],Int[xlist[[m]][[n]][[i]],m,n]},{i,1,Steps}],{n,15,40},{m,1,10}];
Export["\[CapitalDelta]tab.dat",\[CapitalDelta]tab];
Clear[\[CapitalDelta]tab];


\[CapitalDelta]diftab=ParallelTable[Table[{xlist[[m]][[n]][[i]],nsol[[m]][[n]][xlist[[m]][[n]][[i]]]},{i,1,Steps}],{n,15,40},{m,1,10}];


sumerr2=ParallelTable[Sum[Abs[(\[CapitalDelta]tab2[[m]][[n]][[i]][[2]]-\[CapitalDelta]diftab[[m]][[n]][[i]][[2]])/\[CapitalDelta]diftab[[m]][[n]][[i]][[2]]],{i,1,Steps}],{n,1,25},{m,1,10}];
Export["sumerr2.dat",sumerr2];
Clear[sumerr2];


sumerr=ParallelTable[Sum[Abs[(\[CapitalDelta]tab[[m]][[n]][[i]][[2]]-\[CapitalDelta]diftab[[m]][[n]][[i]][[2]])/\[CapitalDelta]diftab[[m]][[n]][[i]][[2]]],{i,1,Steps}],{n,1,25},{m,1,10}];
Export["sumerr.dat",sumerr];
Clear[sumerr];
