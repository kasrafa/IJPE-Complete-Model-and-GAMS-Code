$title Single-scenario, three-period cast-iron closed-loop supply-chain model

$onText
===============================================================================
SCOPE OF THIS PUBLIC REPRODUCIBILITY FILE
===============================================================================
This is a self-contained GAMS implementation of the final corrected model.
It uses exactly one scenario (n1) and three periods (t1, t2, t3). The scenario
index n is intentionally retained on every scenario-dependent flow, inventory,
demand, cost, and feasibility variable. Hence, the code remains traceable to
the mathematical formulation and can later be extended by adding scenarios.

The model retains all three original objectives: economic cost is minimized,
environmental impact is minimized, and employment is maximized. AUGMECON
converts them into a sequence of single-objective MIP subproblems. Economic cost
is the primary objective; environmental and social targets are epsilon
constraints with normalized augmentation terms. A payoff table is generated
first, followed by a two-dimensional epsilon grid and a filtered Pareto set.
Because |N|=1, the mean-deviation risk term is identically zero and is omitted.

All numerical inputs required by this run are included below. Scenario n1,
demand weights, composition coefficients, capacities, costs, return shares,
IoT data, utility data, emission factors, safety stocks, and solver settings
are the reduced one-member instance prepared from the project data. Before any
decision equation is generated, the code validates the fixed inputs and assigns
numerical values to the inventory and activation bounds. Those prepared values
are then used as coefficients in the optimization model.



===============================================================================
UNIT CONVENTION
===============================================================================
Material and product flows: metric tonne (t).
Time: one model planning period.
Money: Australian dollar (AUD), used consistently with the workbook.
Environmental indicator: kg CO2-eq.
Training: worker-hour. Employment: jobs.
$offText

option limrow=0, limcol=0, solprint=off, decimals=6;

*=============================================================================*
* SETS: one scenario, three periods, and one member per physical model set
*=============================================================================*
Sets
    s       raw-material supplier                   / s1 /
    p       candidate iron-ingot plant              / p1 /
    c       candidate cast-iron plant               / c1 /
    d1      automotive demand point                 / d1_1 /
    d2      kitchen-appliance demand point          / d2_1 /
    d3      architecture-decoration demand point    / d3_1 /
    k       existing collection centre              / k1 /
    f       candidate refurbishing centre           / f1 /
    m       candidate remanufacturing centre        / m1 /
    a       candidate recycling centre              / a1 /
    w       candidate disposal centre               / w1 /
    t       planning periods                         / t1*t3 /
    r1      homogeneous iron material               / iron /
    r2      aggregated additive material            / additive /
    u       representative cast-iron product        / gray /
    i       representative IoT alternative          / i1 /
    n       singleton scenario set                  / n1 /
    q       raw scenario fields
            / originalID, demD1, demD2, demD3, trSP, trSC, cIron, cOther /
    route   reverse routes                          / rec, ref, rem, dis /
    costItem economic-cost report items
            / purchase, production, operation, transport, fixed,
              utility, inventory, shortagePenalty, total /
    payoffRun single-objective anchor runs / costOpt, envOpt, socialOpt /
    paretoMetric stored payoff and Pareto fields
            / cost, environmental, social, epsilonEnvironmental,
              epsilonSocial, slackEnvironmental, slackSocial,
              marginD1, marginD2, marginD3, solveTime, modelStatus /
    gEnv environmental epsilon-grid points / ge0*ge10 /
    gSoc social epsilon-grid points / gs0*gs10 /
;

Alias (t,tt);
Alias (gEnv,ggEnv);
Alias (gSoc,ggSoc);

*=============================================================================*
* RAW SCENARIO AND DEMAND DATA
*=============================================================================*
Table scenarioData(n,q) selected published reduced scenario
          originalID  demD1      demD2     demD3      trSP    trSC     cIron    cOther
    n1    207         27199.410  4720.048  1285.1843  21.005  14.1802  100.3281 42.8511
;

Parameter scenarioProbability(n) probability of the singleton scenario
    / n1 1.0 /;

Parameters
    nodeWeightD1(d1)  demand-point share for D1
    nodeWeightD2(d2)  demand-point share for D2
    nodeWeightD3(d3)  demand-point share for D3
    periodWeight(t)   period share for D1 and D2
    d3PeriodWeight(t) period share for D3
;

nodeWeightD1(d1)=1;
nodeWeightD2(d2)=1;
nodeWeightD3(d3)=1;
periodWeight('t1')=0.32; periodWeight('t2')=0.33; periodWeight('t3')=0.35;
d3PeriodWeight('t1')=0.00; d3PeriodWeight('t2')=0.48; d3PeriodWeight('t3')=0.52;

Parameters
    demandD1(d1,t,n) aggregate automotive demand by node and period (t)
    demandD2(d2,t,n) aggregate appliance demand by node and period (t)
    demandD3(d3,t,n) aggregate architecture demand by node and period (t)
;

demandD1(d1,t,n)=scenarioData(n,'demD1')*nodeWeightD1(d1)*periodWeight(t);
demandD2(d2,t,n)=scenarioData(n,'demD2')*nodeWeightD2(d2)*periodWeight(t);
demandD3(d3,t,n)=scenarioData(n,'demD3')*nodeWeightD3(d3)*d3PeriodWeight(t);

*=============================================================================*
* PRODUCT COMPOSITION DATA: one representative product (t/t)
*=============================================================================*
Parameter usageR2(r2,u) aggregated additive share in one tonne of gray iron;
usageR2('additive','gray')=0.06495;

Parameter ironShare(u) iron share that closes each product mass balance;
ironShare(u)=1-sum(r2,usageR2(r2,u));

*=============================================================================*
* SCALARS
*=============================================================================*
Scalars
    nu       facility and IoT investment budget (AUD)          / 2200000 /
    theta    penalty per scenario-level feasibility margin     / 100000 /
    lambda   economic mean-cost weight for the reduced run     / 1.00 /
    phi      net modeled iron make-up coefficient (t per t)    / 0.20 /
    beta     installed-base return rate                        / 0.30 /
    delta    usable recycling-output rate                      / 0.50 /
    gamma    usable remanufacturing-output rate                 / 0.70 /
    omega    usable refurbishment-output rate                  / 0.85 /
    kappaW   water carbon-equivalent factor (kg CO2-eq per L)  / 0.00084 /
    augmentationWeight dimensionless AUGMECON coefficient       / 0.0001 /
    epsilonEnvironmental current environmental epsilon bound   / 1000000000 /
    epsilonSocial current minimum social-performance bound     / 0 /
    rangeEnvironmental environmental payoff range              / 1 /
    rangeSocial social payoff range                            / 1 /
    rangeCost economic payoff range                            / 1 /
    idealEnvironmental best environmental anchor
    nadirEnvironmental worst environmental payoff value
    idealSocial best social anchor
    nadirSocial worst social payoff value
    idealCost best economic anchor
    nadirCost worst economic payoff value
    dominanceTolerance numerical Pareto comparison tolerance    / 0.00001 /
    totalPayoffSolveTime cumulative payoff-table solution time   / 0 /
    totalParetoSolveTime cumulative AUGMECON grid solution time  / 0 /
    numberOfParetoPoints number of unique nondominated points    / 0 /
;

*=============================================================================*
* PURCHASE AND SUPPLIER-TRANSPORT COSTS FOR n1
*=============================================================================*
Parameters
    routeFactorSP(s,p)       supplier-to-ingot route factor
    routeFactorSC(s,c)       supplier-to-cast route factor
    basePurCostR2(r2,n)      n1 additive base purchase cost (AUD per t)
    purCostR1(r1,s,p,n)      iron purchase cost (AUD per t)
    purCostR2(r2,s,c,n)      additive purchase cost (AUD per t)
    trCostR1(r1,s,p,n)       supplier-to-ingot transport coefficient
    trCostR2(r2,s,c,n)       supplier-to-cast transport coefficient
;

routeFactorSP(s,p)=0.90+0.05*(ord(s)-1)+0.05*(ord(p)-1);
routeFactorSC(s,c)=0.90+0.05*(ord(s)-1)+0.05*(ord(c)-1);

* Composition-weighted cost of the five additive categories for gray iron.
* The aggregation preserves the additive purchase cost of the source instance.
basePurCostR2('additive',n)=1347.05079021
    *scenarioData(n,'cOther')/42.8511;

purCostR1(r1,s,p,n)=scenarioData(n,'cIron')*routeFactorSP(s,p);
purCostR2(r2,s,c,n)=basePurCostR2(r2,n)*routeFactorSC(s,c);
trCostR1(r1,s,p,n)=scenarioData(n,'trSP')*routeFactorSP(s,p);
trCostR2(r2,s,c,n)=scenarioData(n,'trSC')*routeFactorSC(s,c);

*=============================================================================*
* CAPACITY DATA (t per period)
*=============================================================================*
Parameters
    capSR1(s,r1) iron supplier capacity
    capSR2(s,r2) additive supplier capacity
    capP(p)      usable ingot-output capacity
    capC(c)      usable fresh-product-output capacity
    capF(f)      usable refurbished-output capacity
    capM(m)      usable remanufactured-output capacity
    capA(a)      usable recycled-iron-output capacity
    capW(w)      disposal-throughput capacity
;

capSR1(s,r1)=5000; capSR2(s,r2)=1000;
capP(p)=7500;
capC(c)=8000;
capF(f)=1800;
capM(m)=1000;
capA(a)=1800;
capW(w)=900;

*=============================================================================*
* PRODUCTION, OPERATION, TRANSPORT, UTILITY, AND INVENTORY COST DATA
*=============================================================================*
Parameters
    prodCostP(p) production cost per t of usable ingot
    prodCostC(c) production cost per t of usable fresh product
    prodCostM(m) production cost per t of usable remanufactured product
    operCostK(k) recovery preparation cost per t sent to rec-ref-rem
    operCostF(f) refurbishing operating cost per t of usable output
    operCostM(m) remanufacturing operating cost per t of usable output
    operCostA(a) recycling operating cost per t of usable output
    operCostW(w) receiving handling and disposal cost per t received
    trCostPC(p,c) cast-ingot transport cost per t
    trCostCD1(c,d1) cast-to-D1 transport cost per t
    trCostCD2(c,d2) cast-to-D2 transport cost per t
    trCostKF(k,f) collection-to-refurbishing transport cost per t
    trCostKM(k,m) collection-to-remanufacturing transport cost per t
    trCostKA(k,a) collection-to-recycling transport cost per t
    trCostKW(k,w) collection-to-disposal transport cost per t
    trCostAP(a,p) recycling-to-ingot transport cost per t
    trCostFD1(f,d1) refurbishing-to-D1 transport cost per t
    trCostFD2(f,d2) refurbishing-to-D2 transport cost per t
    trCostMD3(m,d3) remanufacturing-to-D3 transport cost per t
    utilityCostC(c,t) utility operating cost of an active cast plant
    invCostR1(r1) iron ending-inventory holding cost per t-period
    invCostR2(r2) additive ending-inventory holding cost per t-period
;

prodCostP(p)=125;
prodCostC(c)=180;
prodCostM(m)=145;
operCostK(k)=30;
operCostF(f)=55;
operCostM(m)=70;
operCostA(a)=40;
operCostW(w)=25;

trCostPC(p,c)=8; trCostCD1(c,d1)=16; trCostCD2(c,d2)=13;
trCostKF(k,f)=10; trCostKM(k,m)=12; trCostKA(k,a)=8; trCostKW(k,w)=7;
trCostAP(a,p)=9; trCostFD1(f,d1)=11; trCostFD2(f,d2)=10;
trCostMD3(m,d3)=14;

utilityCostC('c1','t1')=50000; utilityCostC('c1','t2')=51500;
utilityCostC('c1','t3')=53045;

invCostR1(r1)=4; invCostR2(r2)=1.2;

*=============================================================================*
* REVERSE-ALLOCATION DATA
*=============================================================================*
Table muBase(route,k) fraction of each collection stream sent to each route
          k1
    rec   0.30
    ref   0.40
    rem   0.20
    dis   0.10
;

Parameter mu(route,k,t) reverse-allocation fraction;
mu(route,k,t)=muBase(route,k);

*=============================================================================*
* SAFETY STOCK AND TRAINING DATA
*=============================================================================*
Parameters
    safeStockR1(r1,p,t) required ending iron safety stock (t)
    safeStockR2(r2,c,t) required ending additive safety stock (t)
    totalTrainTime(t)   available training time (worker-hour)
;

safeStockR1(r1,p,t)=100;
safeStockR2(r2,c,t)=20;
totalTrainTime(t)=250;

*=============================================================================*
* FACILITY / IoT DATA: fixed cost, jobs, and training time
*=============================================================================*
Parameters
    fixCostP(p,i) fix cost for the ingot plant and IoT option
    fixCostC(c,i) fix cost for the cast plant and IoT option
    fixCostF(f,i) fix cost for the refurbishing centre and IoT option
    fixCostM(m,i) fix cost for the remanufacturing centre and IoT option
    fixCostA(a,i) fix cost for the recycling centre and IoT option
    fixCostW(w,i) fix cost for the disposal centre and IoT option
    jobsP(p,i) jobs created at the ingot plant
    jobsC(c,i) jobs created at the cast plant
    jobsF(f,i) jobs created at the refurbishing centre
    jobsM(m,i) jobs created at the remanufacturing centre
    jobsA(a,i) jobs created at the recycling centre
    jobsW(w,i) jobs created at the disposal centre
    trainP(p,i) training requirement at the ingot plant
    trainC(c,i) training requirement at the cast plant
    trainF(f,i) training requirement at the refurbishing centre
    trainM(m,i) training requirement at the remanufacturing centre
    trainA(a,i) training requirement at the recycling centre
    trainW(w,i) training requirement at the disposal centre
;

fixCostP(p,i)=140000; fixCostC(c,i)=180000;
fixCostF(f,i)=75000;  fixCostM(m,i)=85000;
fixCostA(a,i)=60000;  fixCostW(w,i)=40000;

jobsP(p,i)=15; jobsC(c,i)=18; jobsF(f,i)=8;
jobsM(m,i)=10; jobsA(a,i)=6; jobsW(w,i)=4;

trainP(p,i)=8; trainC(c,i)=10; trainF(f,i)=5;
trainM(m,i)=6; trainA(a,i)=4; trainW(w,i)=3;

*=============================================================================*
* ENVIRONMENTAL DATA USED BY OBJECTIVE O2 AND POST-SOLVE REPORTING
*=============================================================================*
Parameters
    elecConsC(c)  electricity per active cast plant and period (kJ)
    waterConsC(c) water per active cast plant and period (L)
    kappaE(c)     electricity factor (kg CO2-eq per kJ)
    ghgS(s) supplier process-emission factor per t
    ghgP(p) ingot process-emission factor per t usable output
    ghgC(c) cast process-emission factor per t usable output
    ghgF(f) refurbishing process-emission factor per t usable output
    ghgM(m) remanufacturing process-emission factor per t usable output
    ghgA(a) recycling process-emission factor per t usable output
    ghgW(w) disposal process-emission factor per t received
;

elecConsC(c)=25000000;
waterConsC(c)=150000;
kappaE(c)=0.64/3600;
ghgS(s)=0.80;
ghgP(p)=1.50;
ghgC(c)=1.20;
ghgF(f)=0.60;
ghgM(m)=0.70;
ghgA(a)=0.40;
ghgW(w)=0.90;

*=============================================================================*
* INPUT-DATA RULES: checked before any model equation is generated
*=============================================================================*
abort$(card(n)<>1) 'This file must contain exactly one scenario';
abort$(card(t)<>3) 'This file must contain exactly three periods';
abort$(card(s)<>1 or card(p)<>1 or card(c)<>1 or card(d1)<>1
       or card(d2)<>1 or card(d3)<>1 or card(k)<>1 or card(f)<>1
       or card(m)<>1 or card(a)<>1 or card(w)<>1 or card(r2)<>1
       or card(u)<>1 or card(i)<>1)
    'Every non-time, non-scenario index must contain one representative member';
abort$(card(r1)<>1) 'R1 must remain the singleton homogeneous iron set';
abort$(abs(sum(n,scenarioProbability(n))-1)>1e-9)
    'Singleton scenario probability must sum to one';
abort$(smin(n,scenarioProbability(n))<=0)
    'Scenario probabilities must be strictly positive';
abort$(beta<=0 or beta>1 or delta<=0 or delta>1 or gamma<=0 or gamma>1
       or omega<=0 or omega>1 or phi<=0 or phi>1)
    'Return, recovery, and make-up coefficients must lie in (0,1]';
abort$(theta<0 or nu<0)
    'Penalty and investment budget must be nonnegative';
abort$(lambda<0 or lambda>1)
    'The economic mean-cost weight must lie in [0,1]';
abort$(augmentationWeight<=0 or augmentationWeight>=1)
    'The AUGMECON augmentation coefficient must lie in (0,1)';
abort$(abs(sum(d1,nodeWeightD1(d1))-1)>1e-9
       or abs(sum(d2,nodeWeightD2(d2))-1)>1e-9
       or abs(sum(d3,nodeWeightD3(d3))-1)>1e-9
       or abs(sum(t,periodWeight(t))-1)>1e-9
       or abs(sum(t,d3PeriodWeight(t))-1)>1e-9)
    'Demand allocation weights must sum to one';
abort$(smin((d1,t,n),demandD1(d1,t,n))<0
       or smin((d2,t,n),demandD2(d2,t,n))<0
       or smin((d3,t,n),demandD3(d3,t,n))<0)
    'Demand data must be nonnegative';
abort$(smin((n,q),scenarioData(n,q))<0)
    'Raw scenario fields must be nonnegative';
abort$(smin(u,ironShare(u))<=0) 'Every product must have a positive iron share';
abort$(smax(u,abs(ironShare(u)+sum(r2,usageR2(r2,u))-1))>1e-9)
    'Every product composition must sum to one';
abort$(smax((route,k,t),mu(route,k,t))>1 or smin((route,k,t),mu(route,k,t))<0)
    'Reverse-route shares must lie in [0,1]';
abort$(smax((k,t),abs(sum(route,mu(route,k,t))-1))>1e-9)
    'Reverse-route shares must sum to one';
abort$(smax((d3,n),abs(demandD3(d3,'t1',n)))>1e-9)
    'D3 demand must be zero in t1';
abort$(smin((s,r1),capSR1(s,r1))<0 or smin((s,r2),capSR2(s,r2))<0
       or smin(p,capP(p))<0 or smin(c,capC(c))<0
       or smin(f,capF(f))<0 or smin(m,capM(m))<0
       or smin(a,capA(a))<0 or smin(w,capW(w))<0)
    'All supplier, production, recovery, and disposal capacities must be nonnegative';
abort$(smin(p,prodCostP(p))<0 or smin(c,prodCostC(c))<0
       or smin(m,prodCostM(m))<0 or smin(k,operCostK(k))<0
       or smin(f,operCostF(f))<0 or smin(m,operCostM(m))<0
       or smin(a,operCostA(a))<0 or smin(w,operCostW(w))<0
       or smin((c,t),utilityCostC(c,t))<0
       or smin(r1,invCostR1(r1))<0 or smin(r2,invCostR2(r2))<0)
    'Production, operation, utility, and inventory costs must be nonnegative';
abort$(smin((p,c),trCostPC(p,c))<0 or smin((c,d1),trCostCD1(c,d1))<0
       or smin((c,d2),trCostCD2(c,d2))<0 or smin((k,f),trCostKF(k,f))<0
       or smin((k,m),trCostKM(k,m))<0 or smin((k,a),trCostKA(k,a))<0
       or smin((k,w),trCostKW(k,w))<0 or smin((a,p),trCostAP(a,p))<0
       or smin((f,d1),trCostFD1(f,d1))<0 or smin((f,d2),trCostFD2(f,d2))<0
       or smin((m,d3),trCostMD3(m,d3))<0
       or smin((r1,s,p,n),trCostR1(r1,s,p,n))<0
       or smin((r2,s,c,n),trCostR2(r2,s,c,n))<0)
    'All transport coefficients must be nonnegative';
abort$(smin((r1,s,p,n),purCostR1(r1,s,p,n))<0
       or smin((r2,s,c,n),purCostR2(r2,s,c,n))<0)
    'All material purchase costs must be nonnegative';
abort$(smin((r1,p,t),safeStockR1(r1,p,t))<0
       or smin((r2,c,t),safeStockR2(r2,c,t))<0
       or smin(t,totalTrainTime(t))<0)
    'Safety-stock and training-limit inputs must be nonnegative';
abort$(smin((p,i),fixCostP(p,i))<0 or smin((c,i),fixCostC(c,i))<0
       or smin((f,i),fixCostF(f,i))<0 or smin((m,i),fixCostM(m,i))<0
       or smin((a,i),fixCostA(a,i))<0 or smin((w,i),fixCostW(w,i))<0
       or smin((p,i),jobsP(p,i))<0 or smin((c,i),jobsC(c,i))<0
       or smin((f,i),jobsF(f,i))<0 or smin((m,i),jobsM(m,i))<0
       or smin((a,i),jobsA(a,i))<0 or smin((w,i),jobsW(w,i))<0
       or smin((p,i),trainP(p,i))<0 or smin((c,i),trainC(c,i))<0
       or smin((f,i),trainF(f,i))<0 or smin((m,i),trainM(m,i))<0
       or smin((a,i),trainA(a,i))<0 or smin((w,i),trainW(w,i))<0)
    'Facility cost, employment, and training inputs must be nonnegative';
abort$(smin(c,elecConsC(c))<0 or smin(c,waterConsC(c))<0
       or smin(c,kappaE(c))<0 or kappaW<0 or smin(s,ghgS(s))<0
       or smin(p,ghgP(p))<0 or smin(c,ghgC(c))<0 or smin(f,ghgF(f))<0
       or smin(m,ghgM(m))<0 or smin(a,ghgA(a))<0 or smin(w,ghgW(w))<0)
    'Environmental input coefficients must be nonnegative';

Parameters
    supplyCapR1(r1)        total network iron supply per period (B1)
    supplyCapR2(r2)        total network additive supply per period (B1)
    returnPotential(t,n)   total exogenous return potential (B2)
    routeBound(route,t,n)  per-period reverse-route upper bound (B3)
    routeMax(route)        maximum reverse-route bound (B3)
    invUpperR1(r1,p,t)     cumulative iron-inventory upper bound (B4)
    invUpperR2(r2,c,t)     cumulative additive-inventory upper bound (B5)
    bigMP(p) ingot-plant inbound activation bound
    bigMC(c) cast-plant inbound activation bound
;

supplyCapR1(r1)=sum(s,capSR1(s,r1));
supplyCapR2(r2)=sum(s,capSR2(s,r2));

returnPotential(t,n)=0;
returnPotential(t,n)$(ord(t)>1)=beta*(sum(d1,demandD1(d1,t,n))
                                    +sum(d2,demandD2(d2,t,n)));
routeBound(route,t,n)=smax(k,mu(route,k,t))*returnPotential(t,n);
routeMax(route)=smax((t,n),routeBound(route,t,n));

* DERIVED PARAMETERS: these receive numerical values before the model is built.
* Inventory bounds represent cumulative availability, not measured storage size.
invUpperR1(r1,p,t)=ord(t)*supplyCapR1(r1)
    +delta*smax(n,sum(tt$(ord(tt)<=ord(t)),routeBound('rec',tt,n)));
invUpperR2(r2,c,t)=ord(t)*supplyCapR2(r2);

* Finite inbound bounds used by the two forward-facility activation equations.
bigMP(p)=min(sum(r1,supplyCapR1(r1))+delta*routeMax('rec'),
             phi*capP(p)+smax((r1,t),invUpperR1(r1,p,t)));
bigMC(c)=min(sum(p,capP(p))+sum(r2,supplyCapR2(r2)),
             capC(c)+smax(t,sum(r2,invUpperR2(r2,c,t))));

abort$(smax((r1,p,t),safeStockR1(r1,p,t)-invUpperR1(r1,p,t))>1e-9)
    'Iron safety stock exceeds its valid cumulative upper bound';
abort$(smax((r2,c,t),safeStockR2(r2,c,t)-invUpperR2(r2,c,t))>1e-9)
    'Additive safety stock exceeds its valid cumulative upper bound';

*=============================================================================*
* DECISION VARIABLES
*=============================================================================*
Binary Variables
    yP(p,i) ingot-plant and IoT selection
    yC(c,i) cast-plant and IoT selection
    yF(f,i) refurbishing-facility and IoT selection
    yM(m,i) remanufacturing-facility and IoT selection
    yA(a,i) recycling-facility and IoT selection
    yW(w,i) disposal-facility and IoT selection
;

Positive Variables
    xR1SP(r1,s,p,t,n)  iron shipped from suppliers to ingot plants
    xR2SC(r2,s,c,t,n)  additives shipped from suppliers to cast plants
    xPC(p,c,t,n)       usable ingot shipped to cast plants
    xUCD1(u,c,d1,t,n)  fresh product shipped to D1
    xUCD2(u,c,d2,t,n)  fresh product shipped to D2
    xD1K(d1,k,t,n)     D1 returns shipped to collection
    xD2K(d2,k,t,n)     D2 returns shipped to collection
    xKA(k,a,t,n)       collection-to-recycling flow
    xKF(k,f,t,n)       collection-to-refurbishing flow
    xKM(k,m,t,n)       collection-to-remanufacturing flow
    xKW(k,w,t,n)       collection-to-disposal flow
    xR1AP(r1,a,p,t,n)  usable recycled iron shipped to ingot plants
    xUFD1(u,f,d1,t,n)  refurbished product shipped to D1
    xUFD2(u,f,d2,t,n)  refurbished product shipped to D2
    xUMD3(u,m,d3,t,n)  remanufactured product shipped to D3
    invR1(r1,p,t,n)    ending iron inventory
    invR2(r2,c,t,n)    ending additive inventory
    epsD1(n)            shared D1 scenario-level feasibility margin
    epsD2(n)            shared D2 scenario-level feasibility margin
    epsD3(n)            shared D3 scenario-level feasibility margin
    scenarioCost(n)     economic cost before the shortage penalty
    slackEnvironmental  environmental epsilon-constraint slack
    slackSocial         social epsilon-constraint slack
;

Free Variables
    zTotalCost economic objective including shortage penalty
    zEnvironmental environmental objective in kg CO2-eq
    zSocial social objective in jobs
    zAugmented single AUGMECON scalar objective
;

* Period-one reverse-flow initialization is data handling, not a model row.
xD1K.fx(d1,k,'t1',n)=0;
xD2K.fx(d2,k,'t1',n)=0;

*=============================================================================*
* EQUATIONS
*=============================================================================*
Equations
    objectiveDef
    environmentalObjectiveDef
    socialObjectiveDef
    environmentalEpsilonConstr
    socialEpsilonConstr
    augmentedObjectiveDef
    scenarioCostDef(n)
    rawCapR1(r1,s,t,n)
    rawCapR2(r2,s,t,n)
    capPConstr(p,t,n)
    capCConstr(c,t,n)
    capFConstr(f,t,n)
    capMConstr(m,t,n)
    capAConstr(a,t,n)
    capWConstr(w,t,n)
    selectP(p)
    selectC(c)
    selectF(f)
    selectM(m)
    selectA(a)
    selectW(w)
    activateP(p,t,n)
    activateC(c,t,n)
    demandD1First(d1,t,n)
    demandD1Later(d1,t,n)
    demandD2First(d2,t,n)
    demandD2Later(d2,t,n)
    demandD3Later(d3,t,n)
    returnD1Later(d1,t,n)
    returnD2Later(d2,t,n)
    invR1First(r1,p,t,n)
    invR1Later(r1,p,t,n)
    invR2First(r2,c,t,n)
    invR2Later(r2,c,t,n)
    castMassBalance(c,t,n)
    allocRec(k,t,n)
    allocRef(k,t,n)
    allocRem(k,t,n)
    allocDis(k,t,n)
    recyclingYield(a,t,n)
    refurbishingYield(f,t,n)
    remanufacturingYield(m,t,n)
    budgetConstr
    safetyMinR1(r1,p,t,n)
    safetyMaxR1(r1,p,t,n)
    safetyMinR2(r2,c,t,n)
    safetyMaxR2(r2,c,t,n)
    trainingConstr(t)
;

* Equation (1): economic objective. The singleton probability is explicit.
objectiveDef..
    zTotalCost =e= lambda*sum(n,scenarioProbability(n)*scenarioCost(n))
                  +theta*sum(n,scenarioProbability(n)*
                    (epsD1(n)+epsD2(n)+epsD3(n)));

* Equation (4): environmental objective.
* contains exactly the n1 process impacts plus one horizon utility impact.
environmentalObjectiveDef..
    zEnvironmental =e=
      sum((r1,s,p,t,n),ghgS(s)*xR1SP(r1,s,p,t,n))
     +sum((r2,s,c,t,n),ghgS(s)*xR2SC(r2,s,c,t,n))
     +sum((p,c,t,n),ghgP(p)*xPC(p,c,t,n))
     +sum((u,c,d1,t,n),ghgC(c)*xUCD1(u,c,d1,t,n))
     +sum((u,c,d2,t,n),ghgC(c)*xUCD2(u,c,d2,t,n))
     +sum((u,f,d1,t,n),ghgF(f)*xUFD1(u,f,d1,t,n))
     +sum((u,f,d2,t,n),ghgF(f)*xUFD2(u,f,d2,t,n))
     +sum((u,m,d3,t,n),ghgM(m)*xUMD3(u,m,d3,t,n))
     +sum((r1,a,p,t,n),ghgA(a)*xR1AP(r1,a,p,t,n))
     +sum((k,w,t,n),ghgW(w)*xKW(k,w,t,n))
     +sum((c,i,t),yC(c,i)*(kappaE(c)*elecConsC(c)
                            +kappaW*waterConsC(c)));

* Equation (5): social objective measured by jobs at selected facilities.
socialObjectiveDef..
    zSocial =e=
      sum((p,i),jobsP(p,i)*yP(p,i))+sum((c,i),jobsC(c,i)*yC(c,i))
     +sum((f,i),jobsF(f,i)*yF(f,i))+sum((m,i),jobsM(m,i)*yM(m,i))
     +sum((a,i),jobsA(a,i)*yA(a,i))+sum((w,i),jobsW(w,i)*yW(w,i));

* AUGMECON secondary-objective bounds. Positive slacks measure improvement
* beyond each epsilon target and are rewarded only by the tiny normalized term.
environmentalEpsilonConstr..
    zEnvironmental+slackEnvironmental =e= epsilonEnvironmental;
socialEpsilonConstr..
    zSocial-slackSocial =e= epsilonSocial;
augmentedObjectiveDef..
    zAugmented =e= zTotalCost/rangeCost
      -augmentationWeight*(slackEnvironmental/rangeEnvironmental
                           +slackSocial/rangeSocial);

* Complete scenario cost: purchase + production + operation + transport +
* fixed facility/IoT + utility + ending-inventory holding cost. Every
* scenario-dependent transport sum includes the period index t.
scenarioCostDef(n)..
 scenarioCost(n) =e=
    sum((r1,s,p,t),purCostR1(r1,s,p,n)*xR1SP(r1,s,p,t,n))
   +sum((r2,s,c,t),purCostR2(r2,s,c,n)*xR2SC(r2,s,c,t,n))
   +sum((p,c,t),prodCostP(p)*xPC(p,c,t,n))
   +sum((u,c,d1,t),prodCostC(c)*xUCD1(u,c,d1,t,n))
   +sum((u,c,d2,t),prodCostC(c)*xUCD2(u,c,d2,t,n))
   +sum((u,m,d3,t),prodCostM(m)*xUMD3(u,m,d3,t,n))
   +sum((k,a,t),operCostK(k)*xKA(k,a,t,n))
   +sum((k,f,t),operCostK(k)*xKF(k,f,t,n))
   +sum((k,m,t),operCostK(k)*xKM(k,m,t,n))
   +sum((u,f,d1,t),operCostF(f)*xUFD1(u,f,d1,t,n))
   +sum((u,f,d2,t),operCostF(f)*xUFD2(u,f,d2,t,n))
   +sum((u,m,d3,t),operCostM(m)*xUMD3(u,m,d3,t,n))
   +sum((r1,a,p,t),operCostA(a)*xR1AP(r1,a,p,t,n))
   +sum((k,w,t),operCostW(w)*xKW(k,w,t,n))
   +sum((r1,s,p,t),trCostR1(r1,s,p,n)*xR1SP(r1,s,p,t,n))
   +sum((r2,s,c,t),trCostR2(r2,s,c,n)*xR2SC(r2,s,c,t,n))
   +sum((p,c,t),trCostPC(p,c)*xPC(p,c,t,n))
   +sum((u,c,d1,t),trCostCD1(c,d1)*xUCD1(u,c,d1,t,n))
   +sum((u,c,d2,t),trCostCD2(c,d2)*xUCD2(u,c,d2,t,n))
   +sum((k,f,t),trCostKF(k,f)*xKF(k,f,t,n))
   +sum((k,m,t),trCostKM(k,m)*xKM(k,m,t,n))
   +sum((k,a,t),trCostKA(k,a)*xKA(k,a,t,n))
   +sum((k,w,t),trCostKW(k,w)*xKW(k,w,t,n))
   +sum((r1,a,p,t),trCostAP(a,p)*xR1AP(r1,a,p,t,n))
   +sum((u,f,d1,t),trCostFD1(f,d1)*xUFD1(u,f,d1,t,n))
   +sum((u,f,d2,t),trCostFD2(f,d2)*xUFD2(u,f,d2,t,n))
   +sum((u,m,d3,t),trCostMD3(m,d3)*xUMD3(u,m,d3,t,n))
   +sum((p,i),fixCostP(p,i)*yP(p,i))
   +sum((c,i),fixCostC(c,i)*yC(c,i))
   +sum((f,i),fixCostF(f,i)*yF(f,i))
   +sum((m,i),fixCostM(m,i)*yM(m,i))
   +sum((a,i),fixCostA(a,i)*yA(a,i))
   +sum((w,i),fixCostW(w,i)*yW(w,i))
   +sum((c,i,t),utilityCostC(c,t)*yC(c,i))
   +sum((r1,p,t),invCostR1(r1)*invR1(r1,p,t,n))
   +sum((r2,c,t),invCostR2(r2)*invR2(r2,c,t,n));

* Equations (6)-(7): supplier capacities by material, period, and scenario.
rawCapR1(r1,s,t,n)..sum(p,xR1SP(r1,s,p,t,n))=l=capSR1(s,r1);
rawCapR2(r2,s,t,n)..sum(c,xR2SC(r2,s,c,t,n))=l=capSR2(s,r2);

* Equations (8)-(13): usable-output or throughput capacities.
capCConstr(c,t,n)..
    sum((u,d1),xUCD1(u,c,d1,t,n))+sum((u,d2),xUCD2(u,c,d2,t,n))
    =l=capC(c)*sum(i,yC(c,i));
capPConstr(p,t,n)..sum(c,xPC(p,c,t,n))=l=capP(p)*sum(i,yP(p,i));
capFConstr(f,t,n)..
    sum((u,d1),xUFD1(u,f,d1,t,n))+sum((u,d2),xUFD2(u,f,d2,t,n))
    =l=capF(f)*sum(i,yF(f,i));
capMConstr(m,t,n)..sum((u,d3),xUMD3(u,m,d3,t,n))=l=capM(m)*sum(i,yM(m,i));
capAConstr(a,t,n)..sum((r1,p),xR1AP(r1,a,p,t,n))=l=capA(a)*sum(i,yA(a,i));
capWConstr(w,t,n)..sum(k,xKW(k,w,t,n))=l=capW(w)*sum(i,yW(w,i));

* Equation (14), expanded over facility types: at most one IoT alternative.
selectP(p)..sum(i,yP(p,i))=l=1;
selectC(c)..sum(i,yC(c,i))=l=1;
selectF(f)..sum(i,yF(f,i))=l=1;
selectM(m)..sum(i,yM(m,i))=l=1;
selectA(a)..sum(i,yA(a,i))=l=1;
selectW(w)..sum(i,yW(w,i))=l=1;

* Two data-derived inbound-flow activation constraints.
* Recovery-facility activation follows from positive yields and output limits.
* Disposal activation is already imposed by its throughput capacity.
activateP(p,t,n)..
    sum((r1,s),xR1SP(r1,s,p,t,n))+sum((r1,a),xR1AP(r1,a,p,t,n))
    =l=bigMP(p)*sum(i,yP(p,i));
activateC(c,t,n)..
    sum((r2,s),xR2SC(r2,s,c,t,n))+sum(p,xPC(p,c,t,n))
    =l=bigMC(c)*sum(i,yC(c,i));

* Equations (17)-(21): aggregate demand. Epsilon has scenario index n only and is shared
* across all node-period constraints of its demand group.
demandD1First(d1,t,n)$(ord(t)=1)..
    sum((u,c),xUCD1(u,c,d1,t,n))+epsD1(n)=g=demandD1(d1,t,n);
demandD1Later(d1,t,n)$(ord(t)>1)..
    sum((u,c),xUCD1(u,c,d1,t,n))+sum((u,f),xUFD1(u,f,d1,t,n))
    +epsD1(n)=g=demandD1(d1,t,n);
demandD2First(d2,t,n)$(ord(t)=1)..
    sum((u,c),xUCD2(u,c,d2,t,n))+epsD2(n)=g=demandD2(d2,t,n);
demandD2Later(d2,t,n)$(ord(t)>1)..
    sum((u,c),xUCD2(u,c,d2,t,n))+sum((u,f),xUFD2(u,f,d2,t,n))
    +epsD2(n)=g=demandD2(d2,t,n);
demandD3Later(d3,t,n)$(ord(t)>1)..
    sum((u,m),xUMD3(u,m,d3,t,n))+epsD3(n)=g=demandD3(d3,t,n);

* Installed-base return proxy for periods after initialization.
returnD1Later(d1,t,n)$(ord(t)>1)..
    sum(k,xD1K(d1,k,t,n))=e=beta*demandD1(d1,t,n);
returnD2Later(d2,t,n)$(ord(t)>1)..
    sum(k,xD2K(d2,k,t,n))=e=beta*demandD2(d2,t,n);

* Equations (24)-(25): iron balances. phi is multiplied because it is an input
* requirement coefficient. Recycled inflow starts only after t1.
invR1First(r1,p,t,n)$(ord(t)=1)..
    sum(s,xR1SP(r1,s,p,t,n))
    =e=invR1(r1,p,t,n)+phi*sum(c,xPC(p,c,t,n));
invR1Later(r1,p,t,n)$(ord(t)>1)..
    sum(s,xR1SP(r1,s,p,t,n))+sum(a,xR1AP(r1,a,p,t,n))
    +invR1(r1,p,t-1,n)
    =e=invR1(r1,p,t,n)+phi*sum(c,xPC(p,c,t,n));

* Equations (26)-(27): additive balances by additive material.
invR2First(r2,c,t,n)$(ord(t)=1)..
    sum(s,xR2SC(r2,s,c,t,n))
    =e=invR2(r2,c,t,n)
       +sum(u,usageR2(r2,u)*(sum(d1,xUCD1(u,c,d1,t,n))
                             +sum(d2,xUCD2(u,c,d2,t,n))));
invR2Later(r2,c,t,n)$(ord(t)>1)..
    sum(s,xR2SC(r2,s,c,t,n))+invR2(r2,c,t-1,n)
    =e=invR2(r2,c,t,n)
       +sum(u,usageR2(r2,u)*(sum(d1,xUCD1(u,c,d1,t,n))
                             +sum(d2,xUCD2(u,c,d2,t,n))));

* Equation (28): fresh-product iron-content balance, closing total mass with
* the separately tracked additive balances.
castMassBalance(c,t,n)..
    sum(p,xPC(p,c,t,n))
    =e=sum(u,ironShare(u)*(sum(d1,xUCD1(u,c,d1,t,n))
                           +sum(d2,xUCD2(u,c,d2,t,n))));

* Fixed route-allocation shares. Their validated unit sum also implies total
* collection conservation, so no duplicate aggregate balance is generated.
allocRec(k,t,n)..sum(a,xKA(k,a,t,n))
    =e=mu('rec',k,t)*(sum(d1,xD1K(d1,k,t,n))+sum(d2,xD2K(d2,k,t,n)));
allocRef(k,t,n)..sum(f,xKF(k,f,t,n))
    =e=mu('ref',k,t)*(sum(d1,xD1K(d1,k,t,n))+sum(d2,xD2K(d2,k,t,n)));
allocRem(k,t,n)..sum(m,xKM(k,m,t,n))
    =e=mu('rem',k,t)*(sum(d1,xD1K(d1,k,t,n))+sum(d2,xD2K(d2,k,t,n)));
allocDis(k,t,n)..sum(w,xKW(k,w,t,n))
    =e=mu('dis',k,t)*(sum(d1,xD1K(d1,k,t,n))+sum(d2,xD2K(d2,k,t,n)));

* Equations (33)-(35): recovery yields. Unrecovered fractions are process losses.
recyclingYield(a,t,n)..
    sum((r1,p),xR1AP(r1,a,p,t,n))=e=delta*sum(k,xKA(k,a,t,n));
refurbishingYield(f,t,n)..
    sum((u,d1),xUFD1(u,f,d1,t,n))+sum((u,d2),xUFD2(u,f,d2,t,n))
    =e=omega*sum(k,xKF(k,f,t,n));
remanufacturingYield(m,t,n)..
    sum((u,d3),xUMD3(u,m,d3,t,n))=e=gamma*sum(k,xKM(k,m,t,n));

* Equation (36): facility and IoT budget.
budgetConstr..
    sum((p,i),fixCostP(p,i)*yP(p,i))+sum((c,i),fixCostC(c,i)*yC(c,i))
   +sum((f,i),fixCostF(f,i)*yF(f,i))+sum((m,i),fixCostM(m,i)*yM(m,i))
   +sum((a,i),fixCostA(a,i)*yA(a,i))+sum((w,i),fixCostW(w,i)*yW(w,i))
    =l=nu;

* Equations (37)-(38), expanded into lower and upper rows in GAMS.
* A closed plant has zero inventory; an open plant carries at least safety stock.
safetyMinR1(r1,p,t,n)..
    invR1(r1,p,t,n)=g=safeStockR1(r1,p,t)*sum(i,yP(p,i));
safetyMaxR1(r1,p,t,n)..
    invR1(r1,p,t,n)=l=invUpperR1(r1,p,t)*sum(i,yP(p,i));
safetyMinR2(r2,c,t,n)..
    invR2(r2,c,t,n)=g=safeStockR2(r2,c,t)*sum(i,yC(c,i));
safetyMaxR2(r2,c,t,n)..
    invR2(r2,c,t,n)=l=invUpperR2(r2,c,t)*sum(i,yC(c,i));

* Equation (39): workforce-training limit. Selections are shared by periods.
trainingConstr(t)..
    sum((p,i),trainP(p,i)*yP(p,i))+sum((c,i),trainC(c,i)*yC(c,i))
   +sum((f,i),trainF(f,i)*yF(f,i))+sum((m,i),trainM(m,i)*yM(m,i))
   +sum((a,i),trainA(a,i)*yA(a,i))+sum((w,i),trainW(w,i)*yW(w,i))
    =l=totalTrainTime(t);

*=============================================================================*
* PAYOFF TABLE, AUGMECON GRID, AND PARETO SET
*=============================================================================*
Parameters
    payoffTable(payoffRun,paretoMetric) single-objective payoff table
    paretoPoints(gEnv,gSoc,paretoMetric) all optimal AUGMECON grid solutions
    paretoYP(gEnv,gSoc,p,i) facility decisions at each grid solution
    paretoYC(gEnv,gSoc,c,i) facility decisions at each grid solution
    paretoYF(gEnv,gSoc,f,i) facility decisions at each grid solution
    paretoYM(gEnv,gSoc,m,i) facility decisions at each grid solution
    paretoYA(gEnv,gSoc,a,i) facility decisions at each grid solution
    paretoYW(gEnv,gSoc,w,i) facility decisions at each grid solution
    paretoXR1SP(gEnv,gSoc,r1,s,p,t,n) iron-supply flows by grid solution
    paretoXR2SC(gEnv,gSoc,r2,s,c,t,n) additive-supply flows by grid solution
    paretoXPC(gEnv,gSoc,p,c,t,n) ingot flows by grid solution
    paretoXUCD1(gEnv,gSoc,u,c,d1,t,n) fresh D1 flows by grid solution
    paretoXUCD2(gEnv,gSoc,u,c,d2,t,n) fresh D2 flows by grid solution
    paretoXD1K(gEnv,gSoc,d1,k,t,n) D1 return flows by grid solution
    paretoXD2K(gEnv,gSoc,d2,k,t,n) D2 return flows by grid solution
    paretoXKA(gEnv,gSoc,k,a,t,n) recycling-route flows by grid solution
    paretoXKF(gEnv,gSoc,k,f,t,n) refurbishing-route flows by grid solution
    paretoXKM(gEnv,gSoc,k,m,t,n) remanufacturing-route flows by grid solution
    paretoXKW(gEnv,gSoc,k,w,t,n) disposal-route flows by grid solution
    paretoXR1AP(gEnv,gSoc,r1,a,p,t,n) recycled-iron flows by grid solution
    paretoXUFD1(gEnv,gSoc,u,f,d1,t,n) refurbished D1 flows by grid solution
    paretoXUFD2(gEnv,gSoc,u,f,d2,t,n) refurbished D2 flows by grid solution
    paretoXUMD3(gEnv,gSoc,u,m,d3,t,n) remanufactured D3 flows by grid solution
    paretoInvR1(gEnv,gSoc,r1,p,t,n) iron inventory by grid solution
    paretoInvR2(gEnv,gSoc,r2,c,t,n) additive inventory by grid solution
    paretoEpsD1(gEnv,gSoc,n) D1 feasibility margin by grid solution
    paretoEpsD2(gEnv,gSoc,n) D2 feasibility margin by grid solution
    paretoEpsD3(gEnv,gSoc,n) D3 feasibility margin by grid solution
;

Sets
    feasibleGrid(gEnv,gSoc) grid combinations solved to proven optimality
    dominatedGrid(gEnv,gSoc) optimal grid solutions dominated by another point
    efficientGrid(gEnv,gSoc) nondominated grid solutions before deduplication
    duplicateGrid(gEnv,gSoc) repeated nondominated objective vectors
    paretoGrid(gEnv,gSoc) unique nondominated Pareto points
;

Model augmeconModel / all /;
option mip=cplex, optca=0, optcr=0, reslim=600;

* Relax the epsilon constraints while computing the three anchor solutions.
epsilonEnvironmental=1000000000;
epsilonSocial=0;
rangeEnvironmental=1;
rangeSocial=1;

solve augmeconModel using mip minimizing zTotalCost;
abort$(augmeconModel.modelstat<>1)
    'The economic anchor did not reach proven optimality',
    augmeconModel.modelstat, augmeconModel.solvestat;
totalPayoffSolveTime=totalPayoffSolveTime+augmeconModel.resusd;
payoffTable('costOpt','cost')=zTotalCost.l;
payoffTable('costOpt','environmental')=zEnvironmental.l;
payoffTable('costOpt','social')=zSocial.l;
payoffTable('costOpt','marginD1')=sum(n,epsD1.l(n));
payoffTable('costOpt','marginD2')=sum(n,epsD2.l(n));
payoffTable('costOpt','marginD3')=sum(n,epsD3.l(n));
payoffTable('costOpt','solveTime')=augmeconModel.resusd;
payoffTable('costOpt','modelStatus')=augmeconModel.modelstat;

solve augmeconModel using mip minimizing zEnvironmental;
abort$(augmeconModel.modelstat<>1)
    'The environmental anchor did not reach proven optimality',
    augmeconModel.modelstat, augmeconModel.solvestat;
totalPayoffSolveTime=totalPayoffSolveTime+augmeconModel.resusd;
payoffTable('envOpt','cost')=zTotalCost.l;
payoffTable('envOpt','environmental')=zEnvironmental.l;
payoffTable('envOpt','social')=zSocial.l;
payoffTable('envOpt','marginD1')=sum(n,epsD1.l(n));
payoffTable('envOpt','marginD2')=sum(n,epsD2.l(n));
payoffTable('envOpt','marginD3')=sum(n,epsD3.l(n));
payoffTable('envOpt','solveTime')=augmeconModel.resusd;
payoffTable('envOpt','modelStatus')=augmeconModel.modelstat;

solve augmeconModel using mip maximizing zSocial;
abort$(augmeconModel.modelstat<>1)
    'The social anchor did not reach proven optimality',
    augmeconModel.modelstat, augmeconModel.solvestat;
totalPayoffSolveTime=totalPayoffSolveTime+augmeconModel.resusd;
payoffTable('socialOpt','cost')=zTotalCost.l;
payoffTable('socialOpt','environmental')=zEnvironmental.l;
payoffTable('socialOpt','social')=zSocial.l;
payoffTable('socialOpt','marginD1')=sum(n,epsD1.l(n));
payoffTable('socialOpt','marginD2')=sum(n,epsD2.l(n));
payoffTable('socialOpt','marginD3')=sum(n,epsD3.l(n));
payoffTable('socialOpt','solveTime')=augmeconModel.resusd;
payoffTable('socialOpt','modelStatus')=augmeconModel.modelstat;

idealCost=smin(payoffRun,payoffTable(payoffRun,'cost'));
nadirCost=smax(payoffRun,payoffTable(payoffRun,'cost'));
idealEnvironmental=smin(payoffRun,payoffTable(payoffRun,'environmental'));
nadirEnvironmental=smax(payoffRun,payoffTable(payoffRun,'environmental'));
idealSocial=smax(payoffRun,payoffTable(payoffRun,'social'));
nadirSocial=smin(payoffRun,payoffTable(payoffRun,'social'));

rangeEnvironmental=max(nadirEnvironmental-idealEnvironmental,0.000001);
rangeSocial=max(idealSocial-nadirSocial,0.000001);
rangeCost=max(nadirCost-idealCost,0.000001);

* Each grid run is a SINGLE scalar MIP: minimize economic cost while imposing
* epsilon bounds on environmental impact and employment. The normalized slack
* reward removes weakly efficient solutions without changing the primary goal.
loop((gEnv,gSoc),
    epsilonEnvironmental=idealEnvironmental
        +(ord(gEnv)-1)*(nadirEnvironmental-idealEnvironmental)
          /max(card(gEnv)-1,1);
    epsilonSocial=nadirSocial
        +(ord(gSoc)-1)*(idealSocial-nadirSocial)
          /max(card(gSoc)-1,1);

    solve augmeconModel using mip minimizing zAugmented;
    totalParetoSolveTime=totalParetoSolveTime+augmeconModel.resusd;

    paretoPoints(gEnv,gSoc,'epsilonEnvironmental')=epsilonEnvironmental;
    paretoPoints(gEnv,gSoc,'epsilonSocial')=epsilonSocial;
    paretoPoints(gEnv,gSoc,'modelStatus')=augmeconModel.modelstat;
    paretoPoints(gEnv,gSoc,'solveTime')=augmeconModel.resusd;

    if(augmeconModel.modelstat=1,
        feasibleGrid(gEnv,gSoc)=yes;
        paretoPoints(gEnv,gSoc,'cost')=zTotalCost.l;
        paretoPoints(gEnv,gSoc,'environmental')=zEnvironmental.l;
        paretoPoints(gEnv,gSoc,'social')=zSocial.l;
        paretoPoints(gEnv,gSoc,'slackEnvironmental')=slackEnvironmental.l;
        paretoPoints(gEnv,gSoc,'slackSocial')=slackSocial.l;
        paretoPoints(gEnv,gSoc,'marginD1')=sum(n,epsD1.l(n));
        paretoPoints(gEnv,gSoc,'marginD2')=sum(n,epsD2.l(n));
        paretoPoints(gEnv,gSoc,'marginD3')=sum(n,epsD3.l(n));
        paretoYP(gEnv,gSoc,p,i)=yP.l(p,i);
        paretoYC(gEnv,gSoc,c,i)=yC.l(c,i);
        paretoYF(gEnv,gSoc,f,i)=yF.l(f,i);
        paretoYM(gEnv,gSoc,m,i)=yM.l(m,i);
        paretoYA(gEnv,gSoc,a,i)=yA.l(a,i);
        paretoYW(gEnv,gSoc,w,i)=yW.l(w,i);
        paretoXR1SP(gEnv,gSoc,r1,s,p,t,n)=xR1SP.l(r1,s,p,t,n);
        paretoXR2SC(gEnv,gSoc,r2,s,c,t,n)=xR2SC.l(r2,s,c,t,n);
        paretoXPC(gEnv,gSoc,p,c,t,n)=xPC.l(p,c,t,n);
        paretoXUCD1(gEnv,gSoc,u,c,d1,t,n)=xUCD1.l(u,c,d1,t,n);
        paretoXUCD2(gEnv,gSoc,u,c,d2,t,n)=xUCD2.l(u,c,d2,t,n);
        paretoXD1K(gEnv,gSoc,d1,k,t,n)=xD1K.l(d1,k,t,n);
        paretoXD2K(gEnv,gSoc,d2,k,t,n)=xD2K.l(d2,k,t,n);
        paretoXKA(gEnv,gSoc,k,a,t,n)=xKA.l(k,a,t,n);
        paretoXKF(gEnv,gSoc,k,f,t,n)=xKF.l(k,f,t,n);
        paretoXKM(gEnv,gSoc,k,m,t,n)=xKM.l(k,m,t,n);
        paretoXKW(gEnv,gSoc,k,w,t,n)=xKW.l(k,w,t,n);
        paretoXR1AP(gEnv,gSoc,r1,a,p,t,n)=xR1AP.l(r1,a,p,t,n);
        paretoXUFD1(gEnv,gSoc,u,f,d1,t,n)=xUFD1.l(u,f,d1,t,n);
        paretoXUFD2(gEnv,gSoc,u,f,d2,t,n)=xUFD2.l(u,f,d2,t,n);
        paretoXUMD3(gEnv,gSoc,u,m,d3,t,n)=xUMD3.l(u,m,d3,t,n);
        paretoInvR1(gEnv,gSoc,r1,p,t,n)=invR1.l(r1,p,t,n);
        paretoInvR2(gEnv,gSoc,r2,c,t,n)=invR2.l(r2,c,t,n);
        paretoEpsD1(gEnv,gSoc,n)=epsD1.l(n);
        paretoEpsD2(gEnv,gSoc,n)=epsD2.l(n);
        paretoEpsD3(gEnv,gSoc,n)=epsD3.l(n);
    );
);

* Remove dominated points. Cost and environmental impact are minimized; social
* performance is maximized. A second pass removes duplicate objective vectors.
dominatedGrid(gEnv,gSoc)=no;
dominatedGrid(gEnv,gSoc)$feasibleGrid(gEnv,gSoc)=yes$(
    sum((ggEnv,ggSoc)$(
        feasibleGrid(ggEnv,ggSoc)
        and not (sameas(gEnv,ggEnv) and sameas(gSoc,ggSoc))
        and paretoPoints(ggEnv,ggSoc,'cost')
            <= paretoPoints(gEnv,gSoc,'cost')+dominanceTolerance
        and paretoPoints(ggEnv,ggSoc,'environmental')
            <= paretoPoints(gEnv,gSoc,'environmental')+dominanceTolerance
        and paretoPoints(ggEnv,ggSoc,'social')
            >= paretoPoints(gEnv,gSoc,'social')-dominanceTolerance
        and (
            paretoPoints(ggEnv,ggSoc,'cost')
                < paretoPoints(gEnv,gSoc,'cost')-dominanceTolerance
            or paretoPoints(ggEnv,ggSoc,'environmental')
                < paretoPoints(gEnv,gSoc,'environmental')-dominanceTolerance
            or paretoPoints(ggEnv,ggSoc,'social')
                > paretoPoints(gEnv,gSoc,'social')+dominanceTolerance
        )
    ),1)>0
);

efficientGrid(gEnv,gSoc)=feasibleGrid(gEnv,gSoc);
efficientGrid(gEnv,gSoc)$dominatedGrid(gEnv,gSoc)=no;

duplicateGrid(gEnv,gSoc)=no;
duplicateGrid(gEnv,gSoc)$efficientGrid(gEnv,gSoc)=yes$(
    sum((ggEnv,ggSoc)$(
        efficientGrid(ggEnv,ggSoc)
        and (ord(ggEnv)<ord(gEnv)
             or (ord(ggEnv)=ord(gEnv) and ord(ggSoc)<ord(gSoc)))
        and abs(paretoPoints(ggEnv,ggSoc,'cost')
                -paretoPoints(gEnv,gSoc,'cost'))<=dominanceTolerance
        and abs(paretoPoints(ggEnv,ggSoc,'environmental')
                -paretoPoints(gEnv,gSoc,'environmental'))<=dominanceTolerance
        and abs(paretoPoints(ggEnv,ggSoc,'social')
                -paretoPoints(gEnv,gSoc,'social'))<=dominanceTolerance
    ),1)>0
);

paretoGrid(gEnv,gSoc)=efficientGrid(gEnv,gSoc);
paretoGrid(gEnv,gSoc)$duplicateGrid(gEnv,gSoc)=no;
numberOfParetoPoints=sum((gEnv,gSoc)$paretoGrid(gEnv,gSoc),1);

* Restore the economic anchor so the detailed flow report below has a clear,
* reproducible interpretation. Pareto points remain stored in paretoPoints.
epsilonEnvironmental=1000000000;
epsilonSocial=0;
solve augmeconModel using mip minimizing zTotalCost;
abort$(augmeconModel.modelstat<>1)
    'The final economic reporting solve did not reach proven optimality',
    augmeconModel.modelstat, augmeconModel.solvestat;

*=============================================================================*
* POST-SOLVE REPORTING AND LOGIC CHECKS FOR THE ECONOMIC ANCHOR
*=============================================================================*
Scalars
    reportedEnvironmental post-solve environmental indicator (kg CO2-eq)
    reportedSocial        post-solve employment (jobs)
    solveTimeSeconds      solver resource time
    maximumLogicResidual  largest independently calculated logic residual
;

Parameters
    costBreakdown(costItem) economic cost components
    openP(p,i) selected ingot-plant options
    openC(c,i) selected cast-plant options
    openF(f,i) selected refurbishing-facility options
    openM(m,i) selected remanufacturing-facility options
    openA(a,i) selected recycling-facility options
    openW(w,i) selected disposal-facility options
    freshByProduct(u) total fresh output by product
    refurbishedByProduct(u) total refurbished output by product
    remanufacturedByProduct(u) total remanufactured output by product
    demandSlackD1(d1,t,n) residual slack in D1 demand constraints
    demandSlackD2(d2,t,n) residual slack in D2 demand constraints
    demandSlackD3(d3,t,n) residual slack in D3 demand constraints
    ironBalanceResidual(r1,p,t,n) residual in iron balance
    additiveBalanceResidual(r2,c,t,n) residual in additive balance
    returnResidualD1(d1,t,n) residual in D1 return equation
    returnResidualD2(d2,t,n) residual in D2 return equation
;

costBreakdown('purchase')=
    sum((r1,s,p,t,n),purCostR1(r1,s,p,n)*xR1SP.l(r1,s,p,t,n))
   +sum((r2,s,c,t,n),purCostR2(r2,s,c,n)*xR2SC.l(r2,s,c,t,n));
costBreakdown('production')=
    sum((p,c,t,n),prodCostP(p)*xPC.l(p,c,t,n))
   +sum((u,c,d1,t,n),prodCostC(c)*xUCD1.l(u,c,d1,t,n))
   +sum((u,c,d2,t,n),prodCostC(c)*xUCD2.l(u,c,d2,t,n))
   +sum((u,m,d3,t,n),prodCostM(m)*xUMD3.l(u,m,d3,t,n));
costBreakdown('operation')=
    sum((k,a,t,n),operCostK(k)*xKA.l(k,a,t,n))
   +sum((k,f,t,n),operCostK(k)*xKF.l(k,f,t,n))
   +sum((k,m,t,n),operCostK(k)*xKM.l(k,m,t,n))
   +sum((u,f,d1,t,n),operCostF(f)*xUFD1.l(u,f,d1,t,n))
   +sum((u,f,d2,t,n),operCostF(f)*xUFD2.l(u,f,d2,t,n))
   +sum((u,m,d3,t,n),operCostM(m)*xUMD3.l(u,m,d3,t,n))
   +sum((r1,a,p,t,n),operCostA(a)*xR1AP.l(r1,a,p,t,n))
   +sum((k,w,t,n),operCostW(w)*xKW.l(k,w,t,n));
costBreakdown('transport')=
    sum((r1,s,p,t,n),trCostR1(r1,s,p,n)*xR1SP.l(r1,s,p,t,n))
   +sum((r2,s,c,t,n),trCostR2(r2,s,c,n)*xR2SC.l(r2,s,c,t,n))
   +sum((p,c,t,n),trCostPC(p,c)*xPC.l(p,c,t,n))
   +sum((u,c,d1,t,n),trCostCD1(c,d1)*xUCD1.l(u,c,d1,t,n))
   +sum((u,c,d2,t,n),trCostCD2(c,d2)*xUCD2.l(u,c,d2,t,n))
   +sum((k,f,t,n),trCostKF(k,f)*xKF.l(k,f,t,n))
   +sum((k,m,t,n),trCostKM(k,m)*xKM.l(k,m,t,n))
   +sum((k,a,t,n),trCostKA(k,a)*xKA.l(k,a,t,n))
   +sum((k,w,t,n),trCostKW(k,w)*xKW.l(k,w,t,n))
   +sum((r1,a,p,t,n),trCostAP(a,p)*xR1AP.l(r1,a,p,t,n))
   +sum((u,f,d1,t,n),trCostFD1(f,d1)*xUFD1.l(u,f,d1,t,n))
   +sum((u,f,d2,t,n),trCostFD2(f,d2)*xUFD2.l(u,f,d2,t,n))
   +sum((u,m,d3,t,n),trCostMD3(m,d3)*xUMD3.l(u,m,d3,t,n));
costBreakdown('fixed')=
    sum((p,i),fixCostP(p,i)*yP.l(p,i))+sum((c,i),fixCostC(c,i)*yC.l(c,i))
   +sum((f,i),fixCostF(f,i)*yF.l(f,i))+sum((m,i),fixCostM(m,i)*yM.l(m,i))
   +sum((a,i),fixCostA(a,i)*yA.l(a,i))+sum((w,i),fixCostW(w,i)*yW.l(w,i));
costBreakdown('utility')=sum((c,i,t),utilityCostC(c,t)*yC.l(c,i));
costBreakdown('inventory')=
    sum((r1,p,t,n),invCostR1(r1)*invR1.l(r1,p,t,n))
   +sum((r2,c,t,n),invCostR2(r2)*invR2.l(r2,c,t,n));
costBreakdown('shortagePenalty')=
    theta*sum(n,scenarioProbability(n)*(epsD1.l(n)+epsD2.l(n)+epsD3.l(n)));
costBreakdown('total')=
    lambda*sum(costItem$(not sameas(costItem,'total')
                         and not sameas(costItem,'shortagePenalty')),
               costBreakdown(costItem))
   +costBreakdown('shortagePenalty');

reportedEnvironmental=zEnvironmental.l;
reportedSocial=zSocial.l;

solveTimeSeconds=augmeconModel.resusd;

openP(p,i)=yP.l(p,i); openC(c,i)=yC.l(c,i);
openF(f,i)=yF.l(f,i); openM(m,i)=yM.l(m,i);
openA(a,i)=yA.l(a,i); openW(w,i)=yW.l(w,i);

freshByProduct(u)=sum((c,d1,t,n),xUCD1.l(u,c,d1,t,n))
                 +sum((c,d2,t,n),xUCD2.l(u,c,d2,t,n));
refurbishedByProduct(u)=sum((f,d1,t,n),xUFD1.l(u,f,d1,t,n))
                       +sum((f,d2,t,n),xUFD2.l(u,f,d2,t,n));
remanufacturedByProduct(u)=sum((m,d3,t,n),xUMD3.l(u,m,d3,t,n));

demandSlackD1(d1,t,n)=sum((u,c),xUCD1.l(u,c,d1,t,n))
    +sum((u,f)$(ord(t)>1),xUFD1.l(u,f,d1,t,n))+epsD1.l(n)-demandD1(d1,t,n);
demandSlackD2(d2,t,n)=sum((u,c),xUCD2.l(u,c,d2,t,n))
    +sum((u,f)$(ord(t)>1),xUFD2.l(u,f,d2,t,n))+epsD2.l(n)-demandD2(d2,t,n);
demandSlackD3(d3,t,n)$(ord(t)>1)=sum((u,m),xUMD3.l(u,m,d3,t,n))
    +epsD3.l(n)-demandD3(d3,t,n);

returnResidualD1(d1,t,n)$(ord(t)>1)=
    sum(k,xD1K.l(d1,k,t,n))-beta*demandD1(d1,t,n);
returnResidualD2(d2,t,n)$(ord(t)>1)=
    sum(k,xD2K.l(d2,k,t,n))-beta*demandD2(d2,t,n);

ironBalanceResidual(r1,p,t,n)$(ord(t)=1)=
    sum(s,xR1SP.l(r1,s,p,t,n))-invR1.l(r1,p,t,n)
    -phi*sum(c,xPC.l(p,c,t,n));
ironBalanceResidual(r1,p,t,n)$(ord(t)>1)=
    sum(s,xR1SP.l(r1,s,p,t,n))+sum(a,xR1AP.l(r1,a,p,t,n))
    +invR1.l(r1,p,t-1,n)-invR1.l(r1,p,t,n)
    -phi*sum(c,xPC.l(p,c,t,n));

additiveBalanceResidual(r2,c,t,n)$(ord(t)=1)=
    sum(s,xR2SC.l(r2,s,c,t,n))-invR2.l(r2,c,t,n)
    -sum(u,usageR2(r2,u)*(sum(d1,xUCD1.l(u,c,d1,t,n))
                          +sum(d2,xUCD2.l(u,c,d2,t,n))));
additiveBalanceResidual(r2,c,t,n)$(ord(t)>1)=
    sum(s,xR2SC.l(r2,s,c,t,n))+invR2.l(r2,c,t-1,n)-invR2.l(r2,c,t,n)
    -sum(u,usageR2(r2,u)*(sum(d1,xUCD1.l(u,c,d1,t,n))
                          +sum(d2,xUCD2.l(u,c,d2,t,n))));

* Independent post-solve audit of every substantive model relation.
maximumLogicResidual=0;
maximumLogicResidual=max(maximumLogicResidual,
    smax((r1,s,t,n),max(0,sum(p,xR1SP.l(r1,s,p,t,n))-capSR1(s,r1))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((r2,s,t,n),max(0,sum(c,xR2SC.l(r2,s,c,t,n))-capSR2(s,r2))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((c,t,n),max(0,sum((u,d1),xUCD1.l(u,c,d1,t,n))
        +sum((u,d2),xUCD2.l(u,c,d2,t,n))-capC(c)*sum(i,yC.l(c,i)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((p,t,n),max(0,sum(c,xPC.l(p,c,t,n))-capP(p)*sum(i,yP.l(p,i)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((f,t,n),max(0,sum((u,d1),xUFD1.l(u,f,d1,t,n))
        +sum((u,d2),xUFD2.l(u,f,d2,t,n))-capF(f)*sum(i,yF.l(f,i)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((m,t,n),max(0,sum((u,d3),xUMD3.l(u,m,d3,t,n))
        -capM(m)*sum(i,yM.l(m,i)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((a,t,n),max(0,sum((r1,p),xR1AP.l(r1,a,p,t,n))
        -capA(a)*sum(i,yA.l(a,i)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((w,t,n),max(0,sum(k,xKW.l(k,w,t,n))-capW(w)*sum(i,yW.l(w,i)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax(p,max(0,sum(i,yP.l(p,i))-1)));
maximumLogicResidual=max(maximumLogicResidual,
    smax(c,max(0,sum(i,yC.l(c,i))-1)));
maximumLogicResidual=max(maximumLogicResidual,
    smax(f,max(0,sum(i,yF.l(f,i))-1)));
maximumLogicResidual=max(maximumLogicResidual,
    smax(m,max(0,sum(i,yM.l(m,i))-1)));
maximumLogicResidual=max(maximumLogicResidual,
    smax(a,max(0,sum(i,yA.l(a,i))-1)));
maximumLogicResidual=max(maximumLogicResidual,
    smax(w,max(0,sum(i,yW.l(w,i))-1)));
maximumLogicResidual=max(maximumLogicResidual,
    smax((p,t,n),max(0,sum((r1,s),xR1SP.l(r1,s,p,t,n))
        +sum((r1,a),xR1AP.l(r1,a,p,t,n))-bigMP(p)*sum(i,yP.l(p,i)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((c,t,n),max(0,sum((r2,s),xR2SC.l(r2,s,c,t,n))
        +sum(p,xPC.l(p,c,t,n))-bigMC(c)*sum(i,yC.l(c,i)))));
maximumLogicResidual=max(maximumLogicResidual,
    max(0,-smin((d1,t,n),demandSlackD1(d1,t,n))));
maximumLogicResidual=max(maximumLogicResidual,
    max(0,-smin((d2,t,n),demandSlackD2(d2,t,n))));
maximumLogicResidual=max(maximumLogicResidual,
    max(0,-smin((d3,t,n)$(ord(t)>1),demandSlackD3(d3,t,n))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((d1,t,n)$(ord(t)>1),abs(returnResidualD1(d1,t,n))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((d2,t,n)$(ord(t)>1),abs(returnResidualD2(d2,t,n))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((r1,p,t,n),abs(ironBalanceResidual(r1,p,t,n))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((r2,c,t,n),abs(additiveBalanceResidual(r2,c,t,n))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((c,t,n),abs(sum(p,xPC.l(p,c,t,n))
        -sum(u,ironShare(u)*(sum(d1,xUCD1.l(u,c,d1,t,n))
                            +sum(d2,xUCD2.l(u,c,d2,t,n)))))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((k,t,n),abs(sum(a,xKA.l(k,a,t,n))-mu('rec',k,t)
        *(sum(d1,xD1K.l(d1,k,t,n))+sum(d2,xD2K.l(d2,k,t,n))))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((k,t,n),abs(sum(f,xKF.l(k,f,t,n))-mu('ref',k,t)
        *(sum(d1,xD1K.l(d1,k,t,n))+sum(d2,xD2K.l(d2,k,t,n))))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((k,t,n),abs(sum(m,xKM.l(k,m,t,n))-mu('rem',k,t)
        *(sum(d1,xD1K.l(d1,k,t,n))+sum(d2,xD2K.l(d2,k,t,n))))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((k,t,n),abs(sum(w,xKW.l(k,w,t,n))-mu('dis',k,t)
        *(sum(d1,xD1K.l(d1,k,t,n))+sum(d2,xD2K.l(d2,k,t,n))))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((a,t,n),abs(sum((r1,p),xR1AP.l(r1,a,p,t,n))
        -delta*sum(k,xKA.l(k,a,t,n)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((f,t,n),abs(sum((u,d1),xUFD1.l(u,f,d1,t,n))
        +sum((u,d2),xUFD2.l(u,f,d2,t,n))-omega*sum(k,xKF.l(k,f,t,n)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((m,t,n),abs(sum((u,d3),xUMD3.l(u,m,d3,t,n))
        -gamma*sum(k,xKM.l(k,m,t,n)))));
maximumLogicResidual=max(maximumLogicResidual,
    max(0,sum((p,i),fixCostP(p,i)*yP.l(p,i))
        +sum((c,i),fixCostC(c,i)*yC.l(c,i))
        +sum((f,i),fixCostF(f,i)*yF.l(f,i))
        +sum((m,i),fixCostM(m,i)*yM.l(m,i))
        +sum((a,i),fixCostA(a,i)*yA.l(a,i))
        +sum((w,i),fixCostW(w,i)*yW.l(w,i))-nu));
maximumLogicResidual=max(maximumLogicResidual,
    smax((r1,p,t,n),max(0,safeStockR1(r1,p,t)*sum(i,yP.l(p,i))
        -invR1.l(r1,p,t,n))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((r1,p,t,n),max(0,invR1.l(r1,p,t,n)
        -invUpperR1(r1,p,t)*sum(i,yP.l(p,i)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((r2,c,t,n),max(0,safeStockR2(r2,c,t)*sum(i,yC.l(c,i))
        -invR2.l(r2,c,t,n))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((r2,c,t,n),max(0,invR2.l(r2,c,t,n)
        -invUpperR2(r2,c,t)*sum(i,yC.l(c,i)))));
maximumLogicResidual=max(maximumLogicResidual,
    smax(t,max(0,sum((p,i),trainP(p,i)*yP.l(p,i))
        +sum((c,i),trainC(c,i)*yC.l(c,i))
        +sum((f,i),trainF(f,i)*yF.l(f,i))
        +sum((m,i),trainM(m,i)*yM.l(m,i))
        +sum((a,i),trainA(a,i)*yA.l(a,i))
        +sum((w,i),trainW(w,i)*yW.l(w,i))-totalTrainTime(t))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((d1,k,n),abs(xD1K.l(d1,k,'t1',n))));
maximumLogicResidual=max(maximumLogicResidual,
    smax((d2,k,n),abs(xD2K.l(d2,k,'t1',n))));
maximumLogicResidual=max(maximumLogicResidual,
    abs(costBreakdown('total')-zTotalCost.l));

abort$(abs(costBreakdown('total')-zTotalCost.l)>1e-4)
    'Post-solve cost breakdown does not reconcile to the objective';
abort$(smin((d1,t,n),demandSlackD1(d1,t,n)) < -1e-5)
    'A D1 demand constraint is violated';
abort$(smin((d2,t,n),demandSlackD2(d2,t,n)) < -1e-5)
    'A D2 demand constraint is violated';
abort$(smin((d3,t,n)$(ord(t)>1),demandSlackD3(d3,t,n)) < -1e-5)
    'A D3 demand constraint is violated';
abort$(smax((d1,t,n)$(ord(t)>1),abs(returnResidualD1(d1,t,n))) > 1e-5)
    'A D1 return equation is violated';
abort$(smax((d2,t,n)$(ord(t)>1),abs(returnResidualD2(d2,t,n))) > 1e-5)
    'A D2 return equation is violated';
abort$(smax((r1,p,t,n),abs(ironBalanceResidual(r1,p,t,n))) > 1e-5)
    'An iron inventory balance is violated';
abort$(smax((r2,c,t,n),abs(additiveBalanceResidual(r2,c,t,n))) > 1e-5)
    'An additive inventory balance is violated';
abort$(maximumLogicResidual>1e-4)
    'The independent post-solve audit detected a violated model relation',
    maximumLogicResidual;

display payoffTable, paretoPoints, feasibleGrid, dominatedGrid,
        efficientGrid, duplicateGrid, paretoGrid,
        paretoYP, paretoYC, paretoYF, paretoYM, paretoYA, paretoYW,
        paretoXR1SP, paretoXR2SC, paretoXPC, paretoXUCD1, paretoXUCD2,
        paretoXD1K, paretoXD2K, paretoXKA, paretoXKF, paretoXKM, paretoXKW,
        paretoXR1AP, paretoXUFD1, paretoXUFD2, paretoXUMD3,
        paretoInvR1, paretoInvR2, paretoEpsD1, paretoEpsD2, paretoEpsD3,
        idealCost, nadirCost, idealEnvironmental, nadirEnvironmental,
        idealSocial, nadirSocial, rangeCost, rangeEnvironmental, rangeSocial,
        totalPayoffSolveTime, totalParetoSolveTime, numberOfParetoPoints,
        scenarioData, scenarioProbability, demandD1, demandD2, demandD3,
        usageR2, ironShare, returnPotential, routeBound, routeMax,
        invUpperR1, invUpperR2, bigMP, bigMC,
        zTotalCost.l, zEnvironmental.l, zSocial.l, zAugmented.l,
        scenarioCost.l, epsD1.l, epsD2.l, epsD3.l,
        costBreakdown, reportedEnvironmental, reportedSocial,
        solveTimeSeconds, maximumLogicResidual,
        openP, openC, openF, openM, openA, openW,
        freshByProduct, refurbishedByProduct, remanufacturedByProduct,
        demandSlackD1, demandSlackD2, demandSlackD3,
        returnResidualD1, returnResidualD2,
        ironBalanceResidual, additiveBalanceResidual,
        xR1SP.l, xR2SC.l, xPC.l, xUCD1.l, xUCD2.l,
        xD1K.l, xD2K.l, xKA.l, xKF.l, xKM.l, xKW.l,
        xR1AP.l, xUFD1.l, xUFD2.l, xUMD3.l, invR1.l, invR2.l;

File paretoCsv / Single_Scenario_Three_Period_Pareto.csv /;
put paretoCsv;
put 'envGrid,socialGrid,cost,environmental,social,epsilonEnvironmental,'
    'epsilonSocial,slackEnvironmental,slackSocial,marginD1,marginD2,'
    'marginD3,solveTime' /;
loop((gEnv,gSoc)$paretoGrid(gEnv,gSoc),
    put gEnv.tl, ',', gSoc.tl, ',',
        paretoPoints(gEnv,gSoc,'cost'):0:6, ',',
        paretoPoints(gEnv,gSoc,'environmental'):0:6, ',',
        paretoPoints(gEnv,gSoc,'social'):0:6, ',',
        paretoPoints(gEnv,gSoc,'epsilonEnvironmental'):0:6, ',',
        paretoPoints(gEnv,gSoc,'epsilonSocial'):0:6, ',',
        paretoPoints(gEnv,gSoc,'slackEnvironmental'):0:6, ',',
        paretoPoints(gEnv,gSoc,'slackSocial'):0:6, ',',
        paretoPoints(gEnv,gSoc,'marginD1'):0:6, ',',
        paretoPoints(gEnv,gSoc,'marginD2'):0:6, ',',
        paretoPoints(gEnv,gSoc,'marginD3'):0:6, ',',
        paretoPoints(gEnv,gSoc,'solveTime'):0:6 /;
);
putclose;

File paretoDecisionCsv / Single_Scenario_Three_Period_Pareto_Decisions.csv /;
put paretoDecisionCsv;
put 'envGrid,socialGrid,variable,index1,index2,index3,period,scenario,value' /;
loop((gEnv,gSoc,p,i)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',yP,',p.tl,',',i.tl,',,,,',paretoYP(gEnv,gSoc,p,i):0:6 /;);
loop((gEnv,gSoc,c,i)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',yC,',c.tl,',',i.tl,',,,,',paretoYC(gEnv,gSoc,c,i):0:6 /;);
loop((gEnv,gSoc,f,i)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',yF,',f.tl,',',i.tl,',,,,',paretoYF(gEnv,gSoc,f,i):0:6 /;);
loop((gEnv,gSoc,m,i)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',yM,',m.tl,',',i.tl,',,,,',paretoYM(gEnv,gSoc,m,i):0:6 /;);
loop((gEnv,gSoc,a,i)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',yA,',a.tl,',',i.tl,',,,,',paretoYA(gEnv,gSoc,a,i):0:6 /;);
loop((gEnv,gSoc,w,i)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',yW,',w.tl,',',i.tl,',,,,',paretoYW(gEnv,gSoc,w,i):0:6 /;);
loop((gEnv,gSoc,r1,s,p,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xR1SP,',r1.tl,',',s.tl,',',p.tl,',',t.tl,',',n.tl,',',
        paretoXR1SP(gEnv,gSoc,r1,s,p,t,n):0:6 /;);
loop((gEnv,gSoc,r2,s,c,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xR2SC,',r2.tl,',',s.tl,',',c.tl,',',t.tl,',',n.tl,',',
        paretoXR2SC(gEnv,gSoc,r2,s,c,t,n):0:6 /;);
loop((gEnv,gSoc,p,c,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xPC,',p.tl,',',c.tl,',,',t.tl,',',n.tl,',',
        paretoXPC(gEnv,gSoc,p,c,t,n):0:6 /;);
loop((gEnv,gSoc,u,c,d1,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xUCD1,',u.tl,',',c.tl,',',d1.tl,',',t.tl,',',n.tl,',',
        paretoXUCD1(gEnv,gSoc,u,c,d1,t,n):0:6 /;);
loop((gEnv,gSoc,u,c,d2,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xUCD2,',u.tl,',',c.tl,',',d2.tl,',',t.tl,',',n.tl,',',
        paretoXUCD2(gEnv,gSoc,u,c,d2,t,n):0:6 /;);
loop((gEnv,gSoc,d1,k,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xD1K,',d1.tl,',',k.tl,',,',t.tl,',',n.tl,',',
        paretoXD1K(gEnv,gSoc,d1,k,t,n):0:6 /;);
loop((gEnv,gSoc,d2,k,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xD2K,',d2.tl,',',k.tl,',,',t.tl,',',n.tl,',',
        paretoXD2K(gEnv,gSoc,d2,k,t,n):0:6 /;);
loop((gEnv,gSoc,k,a,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xKA,',k.tl,',',a.tl,',,',t.tl,',',n.tl,',',
        paretoXKA(gEnv,gSoc,k,a,t,n):0:6 /;);
loop((gEnv,gSoc,k,f,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xKF,',k.tl,',',f.tl,',,',t.tl,',',n.tl,',',
        paretoXKF(gEnv,gSoc,k,f,t,n):0:6 /;);
loop((gEnv,gSoc,k,m,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xKM,',k.tl,',',m.tl,',,',t.tl,',',n.tl,',',
        paretoXKM(gEnv,gSoc,k,m,t,n):0:6 /;);
loop((gEnv,gSoc,k,w,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xKW,',k.tl,',',w.tl,',,',t.tl,',',n.tl,',',
        paretoXKW(gEnv,gSoc,k,w,t,n):0:6 /;);
loop((gEnv,gSoc,r1,a,p,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xR1AP,',r1.tl,',',a.tl,',',p.tl,',',t.tl,',',n.tl,',',
        paretoXR1AP(gEnv,gSoc,r1,a,p,t,n):0:6 /;);
loop((gEnv,gSoc,u,f,d1,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xUFD1,',u.tl,',',f.tl,',',d1.tl,',',t.tl,',',n.tl,',',
        paretoXUFD1(gEnv,gSoc,u,f,d1,t,n):0:6 /;);
loop((gEnv,gSoc,u,f,d2,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xUFD2,',u.tl,',',f.tl,',',d2.tl,',',t.tl,',',n.tl,',',
        paretoXUFD2(gEnv,gSoc,u,f,d2,t,n):0:6 /;);
loop((gEnv,gSoc,u,m,d3,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',xUMD3,',u.tl,',',m.tl,',',d3.tl,',',t.tl,',',n.tl,',',
        paretoXUMD3(gEnv,gSoc,u,m,d3,t,n):0:6 /;);
loop((gEnv,gSoc,r1,p,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',invR1,',r1.tl,',',p.tl,',,',t.tl,',',n.tl,',',
        paretoInvR1(gEnv,gSoc,r1,p,t,n):0:6 /;);
loop((gEnv,gSoc,r2,c,t,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',invR2,',r2.tl,',',c.tl,',,',t.tl,',',n.tl,',',
        paretoInvR2(gEnv,gSoc,r2,c,t,n):0:6 /;);
loop((gEnv,gSoc,n)$paretoGrid(gEnv,gSoc),
    put gEnv.tl,',',gSoc.tl,',epsD1,',n.tl,',,,,,',paretoEpsD1(gEnv,gSoc,n):0:6 /;
    put gEnv.tl,',',gSoc.tl,',epsD2,',n.tl,',,,,,',paretoEpsD2(gEnv,gSoc,n):0:6 /;
    put gEnv.tl,',',gSoc.tl,',epsD3,',n.tl,',,,,,',paretoEpsD3(gEnv,gSoc,n):0:6 /;
);
putclose;

File payoffCsv / Single_Scenario_Three_Period_Payoff_Table.csv /;
put payoffCsv;
put 'anchor,cost,environmental,social,marginD1,marginD2,marginD3,'
    'solveTime,modelStatus' /;
loop(payoffRun,
    put payoffRun.tl, ',',
        payoffTable(payoffRun,'cost'):0:6, ',',
        payoffTable(payoffRun,'environmental'):0:6, ',',
        payoffTable(payoffRun,'social'):0:6, ',',
        payoffTable(payoffRun,'marginD1'):0:6, ',',
        payoffTable(payoffRun,'marginD2'):0:6, ',',
        payoffTable(payoffRun,'marginD3'):0:6, ',',
        payoffTable(payoffRun,'solveTime'):0:6, ',',
        payoffTable(payoffRun,'modelStatus'):0:0 /;
);
putclose;

File summaryCsv / Single_Scenario_Three_Period_Economic_Anchor.csv /;
put summaryCsv;
put 'metric,value' /;
put 'totalCost,',zTotalCost.l:0:6 /;
put 'environmental,',zEnvironmental.l:0:6 /;
put 'social,',zSocial.l:0:6 /;
put 'marginD1,',sum(n,epsD1.l(n)):0:6 /;
put 'marginD2,',sum(n,epsD2.l(n)):0:6 /;
put 'marginD3,',sum(n,epsD3.l(n)):0:6 /;
put 'solveTimeSeconds,',solveTimeSeconds:0:6 /;
put 'maximumLogicResidual,',maximumLogicResidual:0:6 /;
loop(costItem,
    put 'cost_',costItem.tl,',',costBreakdown(costItem):0:6 /;
);
putclose;

File flowCsv / Single_Scenario_Three_Period_Economic_Anchor_Flows.csv /;
put flowCsv;
put 'variable,index1,index2,index3,period,scenario,value' /;
loop((r1,s,p,t,n)$xR1SP.l(r1,s,p,t,n),
    put 'xR1SP,',r1.tl,',',s.tl,',',p.tl,',',t.tl,',',n.tl,',',
        xR1SP.l(r1,s,p,t,n):0:6 /;);
loop((r2,s,c,t,n)$xR2SC.l(r2,s,c,t,n),
    put 'xR2SC,',r2.tl,',',s.tl,',',c.tl,',',t.tl,',',n.tl,',',
        xR2SC.l(r2,s,c,t,n):0:6 /;);
loop((p,c,t,n)$xPC.l(p,c,t,n),
    put 'xPC,',p.tl,',',c.tl,',,',t.tl,',',n.tl,',',xPC.l(p,c,t,n):0:6 /;);
loop((u,c,d1,t,n)$xUCD1.l(u,c,d1,t,n),
    put 'xUCD1,',u.tl,',',c.tl,',',d1.tl,',',t.tl,',',n.tl,',',
        xUCD1.l(u,c,d1,t,n):0:6 /;);
loop((u,c,d2,t,n)$xUCD2.l(u,c,d2,t,n),
    put 'xUCD2,',u.tl,',',c.tl,',',d2.tl,',',t.tl,',',n.tl,',',
        xUCD2.l(u,c,d2,t,n):0:6 /;);
loop((d1,k,t,n)$xD1K.l(d1,k,t,n),
    put 'xD1K,',d1.tl,',',k.tl,',,',t.tl,',',n.tl,',',xD1K.l(d1,k,t,n):0:6 /;);
loop((d2,k,t,n)$xD2K.l(d2,k,t,n),
    put 'xD2K,',d2.tl,',',k.tl,',,',t.tl,',',n.tl,',',xD2K.l(d2,k,t,n):0:6 /;);
loop((k,a,t,n)$xKA.l(k,a,t,n),
    put 'xKA,',k.tl,',',a.tl,',,',t.tl,',',n.tl,',',xKA.l(k,a,t,n):0:6 /;);
loop((k,f,t,n)$xKF.l(k,f,t,n),
    put 'xKF,',k.tl,',',f.tl,',,',t.tl,',',n.tl,',',xKF.l(k,f,t,n):0:6 /;);
loop((k,m,t,n)$xKM.l(k,m,t,n),
    put 'xKM,',k.tl,',',m.tl,',,',t.tl,',',n.tl,',',xKM.l(k,m,t,n):0:6 /;);
loop((k,w,t,n)$xKW.l(k,w,t,n),
    put 'xKW,',k.tl,',',w.tl,',,',t.tl,',',n.tl,',',xKW.l(k,w,t,n):0:6 /;);
loop((r1,a,p,t,n)$xR1AP.l(r1,a,p,t,n),
    put 'xR1AP,',r1.tl,',',a.tl,',',p.tl,',',t.tl,',',n.tl,',',
        xR1AP.l(r1,a,p,t,n):0:6 /;);
loop((u,f,d1,t,n)$xUFD1.l(u,f,d1,t,n),
    put 'xUFD1,',u.tl,',',f.tl,',',d1.tl,',',t.tl,',',n.tl,',',
        xUFD1.l(u,f,d1,t,n):0:6 /;);
loop((u,f,d2,t,n)$xUFD2.l(u,f,d2,t,n),
    put 'xUFD2,',u.tl,',',f.tl,',',d2.tl,',',t.tl,',',n.tl,',',
        xUFD2.l(u,f,d2,t,n):0:6 /;);
loop((u,m,d3,t,n)$xUMD3.l(u,m,d3,t,n),
    put 'xUMD3,',u.tl,',',m.tl,',',d3.tl,',',t.tl,',',n.tl,',',
        xUMD3.l(u,m,d3,t,n):0:6 /;);
loop((r1,p,t,n)$invR1.l(r1,p,t,n),
    put 'invR1,',r1.tl,',',p.tl,',,',t.tl,',',n.tl,',',invR1.l(r1,p,t,n):0:6 /;);
loop((r2,c,t,n)$invR2.l(r2,c,t,n),
    put 'invR2,',r2.tl,',',c.tl,',,',t.tl,',',n.tl,',',invR2.l(r2,c,t,n):0:6 /;);
putclose;

execute_unload 'Single_Scenario_Three_Period_AUGMECON_Results.gdx',
    scenarioData, scenarioProbability, demandD1, demandD2, demandD3,
    usageR2, ironShare, purCostR1, purCostR2, trCostR1, trCostR2,
    capSR1, capSR2, capP, capC, capF, capM, capA, capW,
    returnPotential, routeBound, routeMax,
    invUpperR1, invUpperR2, bigMP, bigMC,
    yP, yC, yF, yM, yA, yW,
    xR1SP, xR2SC, xPC, xUCD1, xUCD2, xD1K, xD2K,
    xKA, xKF, xKM, xKW, xR1AP, xUFD1, xUFD2, xUMD3,
    invR1, invR2, epsD1, epsD2, epsD3, scenarioCost,
    zTotalCost, zEnvironmental, zSocial, zAugmented,
    costBreakdown, reportedEnvironmental, reportedSocial,
    solveTimeSeconds, maximumLogicResidual,
    payoffTable, paretoPoints, feasibleGrid, dominatedGrid,
    efficientGrid, duplicateGrid, paretoGrid,
    paretoYP, paretoYC, paretoYF, paretoYM, paretoYA, paretoYW,
    paretoXR1SP, paretoXR2SC, paretoXPC, paretoXUCD1, paretoXUCD2,
    paretoXD1K, paretoXD2K, paretoXKA, paretoXKF, paretoXKM, paretoXKW,
    paretoXR1AP, paretoXUFD1, paretoXUFD2, paretoXUMD3,
    paretoInvR1, paretoInvR2, paretoEpsD1, paretoEpsD2, paretoEpsD3,
    idealCost, nadirCost, idealEnvironmental, nadirEnvironmental,
    idealSocial, nadirSocial, rangeCost, rangeEnvironmental, rangeSocial,
    totalPayoffSolveTime, totalParetoSolveTime, numberOfParetoPoints,
    demandSlackD1, demandSlackD2, demandSlackD3,
    returnResidualD1, returnResidualD2,
    ironBalanceResidual, additiveBalanceResidual;
