*cd 'C:/Users/Vianello/Documents/Dropbox/PPIR/PPIR-MA-results'.


*GET 
  FILE='effectsizes.graphdata.sav'.
*DATASET NAME esdata.
*DATASET ACTIVATE esdata.

COMPUTE Overalld=0.
EXECUTE.
if Study='Moral Inversion' Overalld=.51.
if Study='Presumption of Guilt' Overalld=.17.
if Study='Moral Cliff' Overalld=.695.
if Study='Intuitive economics' Overalld=.48.
if Study='Higher Standards - Company' Overalld=.86.
if Study='Higher Standards - Charity' Overalld=.93.
if Study='Cold Hearted Prosociality' Overalld=2.05.
if Study='Burn in Hell' Overalld=.206.
if Study='Bigot-Misanthrope' Overalld=1.27.
if Study='Bad Tipper' Overalld=.568.
if Study='Belief-Act Inconsistency' Overalld=.41.
EXECUTE.


compute LB=0. 
if Study='Moral Inversion' LB=.36.
if Study='Presumption of Guilt' LB=.09.
if Study='Moral Cliff' LB=.637.
if Study='Intuitive economics' LB=.37.
if Study='Higher Standards - Company' LB=.59.
if Study='Higher Standards - Charity' LB=.72.
if Study='Cold Hearted Prosociality' LB=1.86.
if Study='Burn in Hell' LB=.094.
if Study='Bigot-Misanthrope' LB=1.11.
if Study='Bad Tipper' LB=.34.
if Study='Belief-Act Inconsistency' LB=.22.
EXECUTE.

compute UB=0. 
if Study='Moral Inversion' UB=.66.
if Study='Presumption of Guilt' UB=.26.
if Study='Moral Cliff' UB=.752.
if Study='Intuitive economics' UB=.60.
if Study='Higher Standards - Company' UB=1.12.
if Study='Higher Standards - Charity' UB=1.14.
if Study='Cold Hearted Prosociality' UB=2.24.
if Study='Burn in Hell' UB=.318.
if Study='Bigot-Misanthrope' UB=1.42.
if Study='Bad Tipper' UB=.79.
if Study='Belief-Act Inconsistency' UB=.60.
EXECUTE.



COMPUTE OriginalES=0.
EXECUTE.
if Study='Moral Inversion' OriginalES=.79.
if Study='Presumption of Guilt' OriginalES=.03.
if Study='Moral Cliff' OriginalES=.71.
if Study='Intuitive economics' OriginalES=.42.
if Study='Higher Standards - Company' OriginalES=.35.
if Study='Higher Standards - Charity' OriginalES=.93.
if Study='Cold Hearted Prosociality' OriginalES=2.24.
if Study='Burn in Hell' OriginalES=.27.
if Study='Bigot-Misanthrope' OriginalES=.905.
if Study='Bad Tipper' OriginalES=.65.
if Study='Belief-Act Inconsistency' OriginalES=.38.
EXECUTE.

**colored***

*DO IF (UltimoPrincipale=0).
*RECODE LB UB OriginalES Overalld (ELSE=SYSMIS).
*END IF.
*EXECUTE.


GGRAPH 
  /GRAPHDATASET NAME="graphdataset" VARIABLES=Study EffectSize Overalld LB UB OriginalES MISSING= VARIABLEWISE REPORTMISSING=NO 
  /GRAPHSPEC SOURCE=INLINE. 
BEGIN GPL 
SOURCE: s=userSource(id("graphdataset")) 
  DATA: Study=col(source(s), name("Study"), unit.category()) 
  DATA: Overalld=col(source(s), name("Overalld")) 
  DATA: LB=col(source(s), name("LB")) 
  DATA: UB=col(source(s), name("UB")) 
 DATA: EffectSize=col(source(s), name("EffectSize")) 
  DATA: OriginalES=col(source(s), name("OriginalES")) 
  COORD: rect(dim(1,2), transpose()) 
  GUIDE: axis(dim(1))
  GUIDE: axis(dim(2), label("Standardized Mean Difference (d)")) 
 GUIDE: form.line(position(*,0), color(color.black), size(size.".5px"))
GUIDE: legend(aesthetic(aesthetic.color.exterior), label("Sample"))
GUIDE: legend(aesthetic(aesthetic.shape.interior), label("Original Effect Size"))
SCALE:  cat(dim(1), sort.values ("Higher Standards - Charity", "Higher Standards - Company", "Intuitive economics", "Presumption of Guilt", "Bigot-Misanthrope", "Burn in Hell", "Cold Hearted Prosociality", "Moral Cliff",
 "Belief-Act Inconsistency", "Bad Tipper", "Moral Inversion"))
SCALE: linear(dim(2), min(-1.5), max(3))
ELEMENT: point.dodge.symmetric(position(Study*EffectSize), shape.interior(shape.circle), color.interior(color.grey)) 
ELEMENT: interval(position(region.spread.range(Study*(LB+UB))), shape.interior(shape.ibeam), color.interior(color.black), size(size."60%")) 
ELEMENT: point(position((Study*Overalld)), shape.interior(shape.circle), color.interior(color.green), color.exterior(color.black), size(size."2.5%)) 
ELEMENT: point(position((Study*OriginalES)), shape.interior(shape.cross), color.interior(color.magenta), color.exterior(color.lightblue), size(size."3.5%)) 
END GPL.

*ELEMENT: point(position((Study*OriginalES)), shape.interior(shape.square), color.interior(color.blue), color.exterior(color.blue), transparency.interior(transparency."1"), size(size."4%))

*B&W


GGRAPH 
  /GRAPHDATASET NAME="graphdataset" VARIABLES=Study EffectSize Overalld LB UB OriginalES MISSING= VARIABLEWISE REPORTMISSING=NO 
  /GRAPHSPEC SOURCE=INLINE. 
BEGIN GPL 
  SOURCE: s=userSource(id("graphdataset")) 
  DATA: Study=col(source(s), name("Study"), unit.category()) 
  DATA: Overalld=col(source(s), name("Overalld")) 
  DATA: LB=col(source(s), name("LB")) 
  DATA: UB=col(source(s), name("UB")) 
 DATA: EffectSize=col(source(s), name("EffectSize")) 
  DATA: OriginalES=col(source(s), name("OriginalES")) 
  COORD: rect(dim(1,2), transpose()) 
  GUIDE: axis(dim(1)) 
  GUIDE: axis(dim(2), label("Standardized Mean Difference (d)")) 
 GUIDE: form.line(position(*,0), color(color.black), size(size.".5px"))
GUIDE: legend(aesthetic(aesthetic.color.exterior), label("Sample"))
GUIDE: legend(aesthetic(aesthetic.shape.interior), label("Original Effect Size"))
SCALE:  cat(dim(1), sort.values ("Higher Standards - Charity", "Higher Standards - Company", "Intuitive economics", "Presumption of Guilt", "Bigot-Misanthrope", "Burn in Hell", "Cold Hearted Prosociality", "Moral Cliff",
 "Belief-Act Inconsistency", "Bad Tipper", "Moral Inversion"))
SCALE: linear(dim(2), min(-1.5), max(3)) 
ELEMENT: point.dodge.symmetric(position(Study*EffectSize), shape.interior(shape.circle), color.interior(color.grey)) 
ELEMENT: interval(position(region.spread.range(Study*(LB+UB))), shape.interior(shape.ibeam), color.interior(color.black), size(size."60%")) 
ELEMENT: point(position((Study*Overalld)), shape.interior(shape.circle), color.interior(color.black), color.exterior(color.black), size(size."2.5%)) 
ELEMENT: point(position((Study*OriginalES)), shape.interior(shape.circle), color.interior(color.white), color.exterior(color.white), size(size."2%)) 
ELEMENT: point(position((Study*OriginalES)), shape.interior(shape.cross), color.interior(color.gray),  size(size."3.5%)) 
END GPL.



