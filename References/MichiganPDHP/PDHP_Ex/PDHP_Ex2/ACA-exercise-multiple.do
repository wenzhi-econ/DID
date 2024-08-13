*---------------------------------------------------
*      Empirical Exercise
*      Effect of ACA on Health insurance
*---------------------------------------------------
*-----------------------------------------------------------------------------
* You may not need to do this step
* Libraries: bcuse may be use to load some of the dataset from Wooldridge Textbooks
*ssc install bcuse, replace

*Download latest version of csdid and drdid from SSC
*ssc install drdid, replace
*ssc install csdid, replace

* Theme. We will use basic schemes. This step is not necessary

ssc install color_style
set scheme white2

*---------------------------------------------------------------------------------------
*---------------------------------------------------------------------------------------
*-----------------------------------------------------------------------------
* Load data
cd "C:\Users\user\Michigan"
use "PDHP_Ex2\data\ehec_data.dta", clear
*-----------------------------------------------------------------------------
* Do some data manipulations
* Ensure that stfips are numeric
describe
confirm numeric var  stfips

* Make sure year is numeric (and correct)
tab year

* If year of adoption is missing, you need coudl set it to 0
* Remember that missing in Stata is also equivalent to infinity
replace yexp2=0 if yexp2==.

* Create treatment dummy - 1 if treated by that year, 0 otherwise
gen treated     = (year >= yexp2) 
replace treated = 0 if yexp2==0
*-----------------------------------------------------------------------------
* Create variable that identifies a subset of data withour never-treted
gen sample = 1
replace sample =0 if yexp2==0

*---------------------------------------------------------------------------------------
*---------------------------------------------------------------------------------------
*---------------------------------------------------------------------------------------
* Start the analysis
*---------------------------------------------------------------------------------------
* Get TWFE coefficient
* you will need to install reghdfe
ssc install reghdfe
* Multiple fixed effects are added in "abs()"
reghdfe dins treated [aw=W], abs(stfips year)  cluster(stfips)

*---------------------------------------------------------------------------------------
* Get Bacon decomposition 
ssc install bacondecomp, replace
* Stata version allows the use of weights.
bacondecomp dins treated [aw=W], ddetail

* But you can just as well use it without weights
reghdfe dins treated , abs(stfips year)  cluster(stfips)
bacondecomp dins treated , ddetail

*---------------------------------------------------------------------------------------
* Get de Chaisemartin and D'Haultfoeuille Decomposition
ssc install did_multiplegt, replace
* Basic syntax
did_multiplegt dins stfip year treated,  average_effect robust_dynamic  dynamic(2) breps(100) 


*---------------------------------------------------------------------------------------
* Get TWFE event study coefficients

* create event times (rel_year)
gen rel_year = year - yexp
* Set the base at -1
replace rel_year=-1 if yexp==0
* But also set rel_year = -11 to -1
replace rel_year=-1 if rel_year==-11

* We will need a second variable with the values "shifted" +11
gen rel_year2 = rel_year+11
* And those will require to be labeled. Unfortunately Stata doesnt admit
* negative values with factor notation
foreach i of local rel_lab {
	label define rel_year2 `i' "`=`i'-11'", modify
}
levelsof rel_year2, local(rel_lab)
* estimate the TWFE coefficients

* For Plotting, lets use coeffplot, by Ben Jann
ssc install coefplot, replace
coefplot , vertical drop(_cons)

* But could also the matrix produced to create a customized figure
reghdfe dins ib10.rel_year2, abs(stfips year) cluster(stfips) allbase
preserve
	ssc install lbsvmat
	clear
	matrix tbl=r(table)'
	lbsvmat tbl
	gen t=-11+_n
	** dropping constant
	drop in 17
	color_style tableau
	twoway (rbar tbl5 tbl6 t if t<0, pstyle(p1)) (scatter  tbl1 t if t<0, pstyle(p1)) ///
		   (rbar tbl5 tbl6 t if t>=0, pstyle(p2)) (scatter  tbl1 t if t>=0, pstyle(p2)) , ///
		   legend(order(1 "Pre Treatment" 3 "Post Treatment")) ///
		   xtitle("Event time") ytitle("Event Study Estimate") ///
		   note(Note:Static TWFE 0.074(0.013)) yline(0.074)
	graph export PDHP_Ex2\plots\twfe_stata.pdf   , replace 
restore


*---------------------------------------------------------------------------------------
* Callaway and Sant'Anna (2021) procedure
*---------------------------------------------------------------------------------------
* Follow Standard Stata convention. cmd  depvar indepvars [weights] if
csdid dins [w=W], ivar(stfips) ///  name the Panel ID varialbe
				  time(year)   ///  name the time variable
				  gvar(yexp2)  ///  name of the first treatment period variable 
				  method(reg)  ///  estimation method. "drimp" is the default method
				  saverif(rif1) replace // Note that CSDID estimates Asymptotic SE by default. To get Wbootstrap you should "save" the RIFs
** Default is to use Never treated. Otherwise Indicate "notyet"				  
preserve
	use rif1, clear				
	* Default uses 999 repetitions  
	*csdid_stats event,  wboot 
	* to select a window
	csdid_stats event,  wboot  window(-5 5) estore(never)
	
restore

*** Control groups: Not yet treated

* Use not-yet-treated as comparison group
csdid dins [w=W], ivar(stfips) ///  name the Panel ID varialbe
				  time(year)   ///  name the time variable
				  gvar(yexp2)  ///  name of the first treatment period variable 
				  method(reg)  ///  estimation method. "drimp" is the default method
				  notyet      /// Request using NOTyet treated as well as Never treated as controls
				  saverif(rif1) replace // Note that CSDID estimates Asymptotic SE by default. To get Wbootstrap you should "save" the RIFs
** Default is to use Never treated. Otherwise Indicate "notyet"				  
preserve
	use rif1, clear				
	* Default uses 999 repetitions  
	*csdid_stats event,  wboot 
	* to select a window
	csdid_stats event,  wboot  window(-5 5) estore(notyet)
restore

*---------------------------------------------------------------------------------------

* Use not-yet-treated as comparison group (drop all never-treated)
csdid dins [w=W] if yexp2>0, ivar(stfips) ///  name the Panel ID varialbe, excludes the never treated
				  time(year)   ///  name the time variable
				  gvar(yexp2)  ///  name of the first treatment period variable 
				  method(reg)  ///  estimation method. "drimp" is the default method
				  notyet      /// Request using NOTyet treated as well as Never treated as controls
				  saverif(rif1) replace // Note that CSDID estimates Asymptotic SE by default. To get Wbootstrap you should "save" the RIFs
** Default is to use Never treated. Otherwise Indicate "notyet"				  
preserve
	use rif1, clear				
	* Default uses 999 repetitions  
	*csdid_stats event,  wboot 
	* to select a window
	csdid_stats event,  wboot  window(-5 5) estore(nyever)
restore

* Notice that Uniform Confidence Intervals are stored in e(cband)
* Now we are ready to go! Let me put all the outputs into a table

*** Making graphs using Not yet

csdid dins [w=W], ivar(stfips) ///  name the Panel ID varialbe
				  time(year)   ///  name the time variable
				  gvar(yexp2)  ///  name of the first treatment period variable 
				  method(reg)  ///  estimation method. "drimp" is the default method
				  notyet      // Request using NOTyet treated as well as Never treated as controls

** Default is to use Never treated. Otherwise Indicate "notyet"				  
csdid_plot, group(2014) title(Group 2014) name(m14,replace) legend( row(1))
csdid_plot, group(2015) title(Group 2015) name(m15,replace) legend( row(1))
csdid_plot, group(2016) title(Group 2016) name(m16,replace) legend( row(1))
csdid_plot, group(2017) title(Group 2017) name(m17,replace) legend( row(1))
csdid_plot, group(2019) title(Group 2019) name(m19,replace) legend( row(1))

** One legend graph
net install  grc1leg, from(http://www.stata.com/users/vwiggins)
grc1leg m14 m15 m16 m17 m19 , nocopies ycommon
grc1leg m14 m15 m16 m17 m19 , nocopies ycommon xcommon note(note: ATT(g,t) within cohorts)

graph export PDHP_Ex2\plots\cs_attgt.pdf   , replace 

** For other groups

csdid dins [w=W], ivar(stfips) ///  name the Panel ID varialbe
				  time(year)   ///  name the time variable
				  gvar(yexp2)  ///  name of the first treatment period variable 
				  method(reg)  ///  estimation method. "drimp" is the default method
				  saverif(rif1) replace // Note that CSDID estimates Asymptotic SE by default. To get Wbootstrap you should "save" the RIFs
				  
preserve
	use rif1, clear				
	* Default uses 999 repetitions  
	*csdid_stats event,  wboot 
	* to select a window
	csdid_stats event,  wboot  window(-5 5) 
	csdid_plot, name(never) title(Never Treated) 
	
restore		  


csdid dins [w=W]  , ivar(stfips) ///  name the Panel ID varialbe, excludes the never treated
				  time(year)   ///  name the time variable
				  gvar(yexp2)  ///  name of the first treatment period variable 
				  method(reg)  ///  estimation method. "drimp" is the default method
				  notyet      /// Request using NOTyet treated as well as Never treated as controls
				  saverif(rif1) replace // Note that CSDID estimates Asymptotic SE by default. To get Wbootstrap you should "save" the RIFs
				  
preserve
	use rif1, clear				
	* Default uses 999 repetitions  
	*csdid_stats event,  wboot 
	* to select a window
	csdid_stats event,  wboot  window(-5 5) 
	csdid_plot, name(notyet) title(Not yet treated)
restore		  

csdid dins [w=W] if yexp2>0, ivar(stfips) ///  name the Panel ID varialbe, excludes the never treated
				  time(year)   ///  name the time variable
				  gvar(yexp2)  ///  name of the first treatment period variable 
				  method(reg)  ///  estimation method. "drimp" is the default method
				  notyet      /// Request using NOTyet treated as well as Never treated as controls
				  saverif(rif1) replace // Note that CSDID estimates Asymptotic SE by default. To get Wbootstrap you should "save" the RIFs
				  
preserve
	use rif1, clear				
	* Default uses 999 repetitions  
	*csdid_stats event,  wboot 
	* to select a window
	csdid_stats event,  wboot  window(-5 5) 
	csdid_plot, name(notyetever) title(Ever Treated)
restore		  
grc1leg never notyet notyetever , nocopies ycommon col(3)  name(cmbny)
graph combine cmbny, xsize(10)
graph export PDHP_Ex2\plots\cs_nynv.pdf   , replace 

