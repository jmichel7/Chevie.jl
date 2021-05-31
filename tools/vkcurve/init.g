#############################################################################
##
#A  init.g       VKCURVE package         David Bessis and Jean Michel
##
#Y  Copyright (C) 2001 - 2002  University Paris VII, France.
##
##  This is the init file of the VKCURVE package.
## 
#############################################################################

if not IsBound(VKCURVE) then 
VKCURVE:=rec(path:=LOADED_PACKAGES.vkcurve);
fi;

VKCURVE.name:="vkcurve";
VKCURVE.version:="1.2";
VKCURVE.date:=[2009,3];
VKCURVE.homepage:="http://webusers.imj-prg.fr/~jean.michel/vkcurve.html";
VKCURVE.copyright:=
"(C) David Bessis, Jean Michel -- compute Pi_1 of hypersurface complements";
VKCURVE.monodromyApprox:=false;
VKCURVE.showSingularProj:=false; 
VKCURVE.showBraiding:=false;
VKCURVE.showLoops:=false;
VKCURVE.showAction:=false;
VKCURVE.showSegments:=false;
VKCURVE.showInsideSegments:=false;
VKCURVE.showWorst:=false;
VKCURVE.showZeros:=false;
VKCURVE.showNewton:=false;
VKCURVE.showgetbraid:=false;
VKCURVE.showRoots:=false;
VKCURVE.showallnewton:=false;# for NewtonRoot
VKCURVE.NewtonLim:=800;      # for NewtonRoot
VKCURVE.AdaptivityFactor:=10; # for ApproxFollowMonodromy
VKCURVE.shrinkBraid:=false;
VKCURVE.mvp:=1;

PrintPkgInit(VKCURVE);

ReadVK:= function(name)
  if not ReadPath(VKCURVE.path, name, ".g", "ReadVK") then
     Error("VKCURVE library file '", name, "' must exist and be readable");
  fi;
end;

AUTO(ReadVK("action"),VKQuotient,DBVKQuotient,BnActsOnFn);
AUTO(ReadVK("loops"),LoopsAroundPunctures);
AUTO(ReadVK("plbraid"),LBraidToWord);
AUTO(ReadVK("polyroot"),NewtonRoot,SeparateRoots,SeparateRootsInitialGuess,
  FindRoots);
AUTO(ReadVK("pres"),ShrinkPresentation,DisplayPresentation,
 TryConjugatePresentation,ConjugatePresentation);
AUTO(ReadVK("global"),PrepareFundamentalGroup,FinishFundamentalGroup,
  FundamentalGroup);
AUTO(ReadVK("util"), BigNorm, ComplexRational,
  DecimalLog, Discy, Dispersal, DistSeg, ResultantMat, Rho, SmallNorm);
AUTO(ReadVK("truemono"),FollowMonodromy);
AUTO(ReadVK("segtobrd"),ApproxFollowMonodromy);
AUTO(ReadVK("mvp"),MvpOps,Mvp,IsMvp,ScalMvp,FactorQuadraticForm,
Jacobian,Hessian,MonomialGcd,OnPolynomials);
AUTO(ReadVK("mvrf"),RatFrac,IsRatFrac,MvpGcd,MvpLcm,LaurentDenominator);

Horner:=function(arg)
  Print("#W Horner is in the GAP library under the name ValuePol\n");
  return ApplyFunc(ValuePol,arg);
end;

Deriv:=function(arg)
  Print("#W Deriv is obsolete, use Derivative\n");
  return ApplyFunc(Derivative,arg);
end;

PresentationOps.Display:=DisplayPresentation; # so as to be able to do this
ReadVK("mvp");
