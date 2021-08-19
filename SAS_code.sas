
%macro simdata(nsim=, n=, nT=, b0=, b1=, b2=, bC=, c0=, c1=, sigma=, pC=, tauC=, tau0=, tau1=, tau2=, tau3=, tauCA1=, tauCA2=, a0=, a1=, a2=, aC=);
data SEM;
call streaminit(123);       /* set random number seed */
DO sim=1 to &nsim by 1;
	DO ID=1 to &n by 1;
		b0=&b0;
		b1=&b1;
		b2=&b2;
		bC=&bC;
		c0=&c0;
		c1=&c1;
		pC=&pC;
		sigma=&sigma;
		tau0=&tau0;
		tau1=&tau1;
		tau2=&tau2;
		tauC=&tauC;

		if(&tau0>0) then	u0= rand('Normal',0,&tau0);
		else u0=0;

		u1= rand('Normal',0,&tau1);
		u2= rand('Normal',0,&tau2);
		uD0= rand('Normal',0,&tau3);
		
		if(&tauC>0) then	uC= rand('Normal',0,&tauC);
		else uC=0;

		if(&tauCA1>0) then	uMcY= rand('Normal',0,&tauCA1);
		else uMcy=0;

		if(&tauCA2>0) then     uMcA= rand('Normal',0,&tauCA2);
		else uMca=0;

		Y_lag_0=0;
		A_lag_0=0;
		D_lag_0=0;
		Y_lag_1=0;
		A_lag_1=0;
		Y_lag_2=0;
		A_lag_2=0;
		C_lag_1=0;
		C_lag_2=0;

		DO PER=1 to &nT;
		if(PER=1) then 
			do;
			e_0=rand('Normal',0,&sigma);
			Y_lag_0=&b0 + u0 + e_0;

	
		C_lag_0 = rand('Bernoulli',pC);
			if(C_lag_0=0) then do;
			C_lag_0=-&pC;
			end;
			else do;
			C_lag_0=(1-&pC);
			end;
			
			pA = exp(&a0+&aC*(C_lag_0)+(&a1+uMcY)*(Y_lag_0))/(1+exp(&a0+&aC*(C_lag_0)+(&a1+uMcY)*(Y_lag_0)));
			uA = rand('UNIFORM');
			A_lag_0= (pA>uA);
			pD = exp(&c0+uD0+&c1*(Y_lag_0))/(1+exp(&c0+uD0+&c1*(Y_lag_0)));
			uD = rand('UNIFORM');
			D_lag_0= (pD>uD);
			Z_lag_0=exp(Y_lag_0);
			end;

		else if(PER=2) then 
			do;
			C_lag_2=C_lag_1;
			C_lag_1=C_lag_0;
			Y_lag_1=Y_lag_0;
			A_lag_1=A_lag_0;
			e_1=rand('Normal',0,&sigma);

			

			Y_lag_0= &b0 + u0 + (&b1+u1)*A_lag_1 + (&bC+uC)*C_lag_1+ e_1;

			C_lag_0 = rand('Bernoulli',pC);
			if(C_lag_0=0) then do;
			C_lag_0=-&pC;
			end;
			else do;
			C_lag_0=(1-&pC);
			end;

			pA = exp(&a0+&aC*(C_lag_0)+(&a1+uMcY)*(Y_lag_0)+(&a2+uMcA)*A_lag_1)/(1+exp(&a0+&aC*(C_lag_0)+(&a1+uMcY)*(Y_lag_0)+(&a2+uMcA)*A_lag_1));
			uA = rand('UNIFORM');
			A_lag_0= (pA>uA);
			pD = exp(&c0+uD0+&c1*(Y_lag_0))/(1+exp(&c0+uD0+&c1*(Y_lag_0)));
			uD = rand('UNIFORM');
			D_lag_0= (pD>uD);
			Z_lag_0=exp(Y_lag_0);
			end;
		else do;
			C_lag_2=C_lag_1;
			Y_lag_2=Y_lag_1;
			A_lag_2=A_lag_1;
			C_lag_1=C_lag_0;
			Y_lag_1=Y_lag_0;
			A_lag_1=A_lag_0;
			e2 = rand('Normal',0,&sigma);

			Y_lag_0= &b0 + u0 + (&b1+u1)*A_lag_1 + (&b2+u2)*A_lag_2 + (&bC+uC)*C_lag_1+e2;

			C_lag_0 = rand('Bernoulli',pC);
			if(C_lag_0=0) then do;
			C_lag_0=-&pC;
			end;
			else do;
			C_lag_0=(1-&pC);
			end;

			pA = exp(&a0+&aC*(C_lag_0)+(&a1+uMcY)*(Y_lag_0)+(&a2+uMcA)*A_lag_1)/(1+exp(&a0+&aC*(C_lag_0)+(&a1+uMcY)*(Y_lag_0)+(&a2+uMcA)*A_lag_1));
			uA = rand('UNIFORM');
			A_lag_0= (pA>uA);
			pD = exp(&c0+uD0+&c1*(Y_lag_0))/(1+exp(&c0+uD0+&c1*(Y_lag_0)));
			uD = rand('UNIFORM');
			D_lag_0= (pD>uD);
			Z_lag_0=exp(Y_lag_0);
			end;
		OUTPUT;
		END;
	END;	
END;  
run;
%mend simdata;

*EXAMPLE 1;
%simdata(nsim=1000, n=1000, nT=100, b0=120, b1=-10, b2=-5, bC=5, c0=-10, c1=0.10, sigma=1, pC=0.3, tauC=0, tau0=5 , tau1=10 , tau2=5, tau3=1, tauCA1=0, tauCA2=0, a0=-3, a1=0.05, a2=1,aC=0.7); 


*Time varying confounding - issues standard analyses;
ODS HTML CLOSE;
ods output SolutionF=testF_wrong(keep = sim Effect Estimate);
proc mixed data=SEM method=ML ASYCOV;
	class PER;
	where PER<=100;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1 /solution;
	by sim;
run;

ODS HTML;
proc means data=testF_wrong;
where effect='A_lag_1';
var estimate;
run;

proc means data=testF_wrong;
where effect='A_lag_2';
var estimate;
run;

ods output SolutionF=testF_wrong2(keep = sim Effect Estimate);
proc mixed data=SEM method=ML ASYCOV;
	where PER=3;
	class PER;
	model Y_lag_0 = Y_lag_1 A_lag_2 A_lag_1 C_lag_1 /solution;
	by sim;
run;
ODS HTML;

ODS HTML CLOSE;
proc sort data=Sem;
by sim;
run;

data Ex1;
set SEM;
where sim=1;
run;

proc export data = Ex1
outfile ="D:\Documents\PhD\Infering_idividual_causal\SIM.csv"
dbms = csv
replace ;
run;


ods output SolutionF=testF(keep = sim Effect Estimate) CovParms=testCov(keep = sim CovParm Estimate) SolutionR=testR(keep= sim Effect Subject Estimate  rename=(Effect=EffectR Estimate=EstimateR)) G=G V=V;
proc mixed data=SEM method=REML ASYCOV;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB outp=est_eps;;
	random int A_lag_2 A_lag_1 /subject=ID type=VC solution G V;
	by sim;
run;

ods output SolutionF=testFF(keep = sim Effect Estimate) CovParms=testCovF(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where ID<=500;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFC(keep = sim Effect Estimate) CovParms=testCovC(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where ID<=100;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFG(keep = sim Effect Estimate) CovParms=testCovG(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=3;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFD(keep = sim Effect Estimate) CovParms=testCovD(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=3 AND ID<=500;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFA(keep = sim Effect Estimate) CovParms=testCovA(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=3 AND ID<=100;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFH(keep = sim Effect Estimate) CovParms=testCovH(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=10;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFE(keep = sim Effect Estimate) CovParms=testCovE(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=10 AND ID<=500;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFB(keep = sim Effect Estimate) CovParms=testCovB(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=10 AND ID<=100;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

DATA Testcova;
set Testcova;
scenario='A';
run;

DATA Testcovb;
set Testcovb;
scenario='B';
run;

DATA Testcovc;
set Testcovc;
scenario='C';
run;

DATA Testcovd;
set Testcovd;
scenario='D';
run;

DATA Testcove;
set Testcove;
scenario='E';
run;

DATA Testcovf;
set Testcovf;
scenario='F';
run;

DATA Testcovg;
set Testcovg;
scenario='G';
run;

DATA TestcovH;
set TestcovH;
scenario='H';
run;

DATA TestcovI;
set Testcov;
scenario='I';
run;

DATA Testcov;
set Testcova Testcovb Testcovc Testcovd Testcove Testcovf Testcovg Testcovh Testcovi;
run;

proc export data = Testcov
outfile ="D:\Documents\PhD\Infering_idividual_causal\EX1_cov.csv"
dbms = csv
replace ;
run;

*fixed;
DATA Testfa;
set Testfa;
scenario='A';
run;

DATA Testfb;
set Testfb;
scenario='B';
run;

DATA Testfc;
set Testfc;
scenario='C';
run;

DATA Testfd;
set Testfd;
scenario='D';
run;

DATA Testfe;
set Testfe;
scenario='E';
run;

DATA Testff;
set Testff;
scenario='F';
run;

DATA Testfg;
set Testfg;
scenario='G';
run;

DATA TestfH;
set TestfH;
scenario='H';
run;

DATA TestfI;
set Testf;
scenario='I';
run;

DATA Testf;
set Testfa Testfb Testfc Testfd Testfe Testff Testfg Testfh Testfi;
run;

proc export data = Testf
outfile ="D:\Documents\PhD\Infering_idividual_causal\EX1_fixed.csv"
dbms = csv
replace ;
run;

*EXAMPLE 3 (same data as EXAMPLE 1);

proc sort data=Sem;
by sim;
run;

ODS HTML CLOSE;
ods output ParameterEstimates=glimF(keep = sim Effect Estimate) FitStatistics=testFitD CovParms=glimCov(keep = sim CovParm Estimate) SolutionR=testRD(keep= sim Effect Subject Estimate  rename=(Effect=EffectR Estimate=EstimateR));
proc glimmix data=Sem method=LAPLACE;
	class D_lag_0;
	model D_lag_0(event='1') = Y_lag_0 /dist=binary link=logit solution;
	random int/ subject=ID solution;
	output out=est_prob(keep=sim id per pD_est D_lag_0) pred(ILINK)=pD_est;
	by sim;
run;

ods output ParameterEstimates=glimFA(keep = sim Effect Estimate) CovParms=glimCovA(keep = sim CovParm Estimate);
proc glimmix data=Sem method=LAPLACE;
	where PER<=3 AND ID<=100;
	class D_lag_0;
	model D_lag_0(event='1') = Y_lag_0 /dist=binary link=logit solution;
	random int/ subject=ID;
	by sim;
run;

ods output ParameterEstimates=glimFB(keep = sim Effect Estimate) CovParms=glimCovB(keep = sim CovParm Estimate);
proc glimmix data=Sem method=LAPLACE;
	where PER<=3 AND ID<=500;
	class D_lag_0;
	model D_lag_0(event='1') = Y_lag_0 /dist=binary link=logit solution;
	random int/ subject=ID;
	by sim;
run;

ods output ParameterEstimates=glimFC(keep = sim Effect Estimate) CovParms=glimCovC(keep = sim CovParm Estimate);
proc glimmix data=Sem method=LAPLACE;
	where PER<=3;
	class D_lag_0;
	model D_lag_0(event='1') = Y_lag_0 /dist=binary link=logit solution;
	random int/ subject=ID;
	by sim;
run;

ods output ParameterEstimates=glimFD(keep = sim Effect Estimate) CovParms=glimCovD(keep = sim CovParm Estimate);
proc glimmix data=Sem method=LAPLACE;
	where PER<=10 AND ID<=100;
	class D_lag_0;
	model D_lag_0(event='1') = Y_lag_0 /dist=binary link=logit solution;
	random int/ subject=ID;
	by sim;
run;

ods output ParameterEstimates=glimFE(keep = sim Effect Estimate) CovParms=glimCovE(keep = sim CovParm Estimate);
proc glimmix data=Sem method=LAPLACE;
	where PER<=10 AND ID<=500;
	class D_lag_0;
	model D_lag_0(event='1') = Y_lag_0 /dist=binary link=logit solution;
	random int/ subject=ID;
	by sim;
run;

ods output ParameterEstimates=glimFF(keep = sim Effect Estimate) CovParms=glimCovF(keep = sim CovParm Estimate);
proc glimmix data=Sem method=LAPLACE;
	where PER<=10;
	class D_lag_0;
	model D_lag_0(event='1') = Y_lag_0 /dist=binary link=logit solution;
	random int/ subject=ID;
	by sim;
run;

ods output ParameterEstimates=glimFG(keep = sim Effect Estimate) CovParms=glimCovG(keep = sim CovParm Estimate);
proc glimmix data=Sem method=LAPLACE;
	where ID<=100;
	class D_lag_0;
	model D_lag_0(event='1') = Y_lag_0 /dist=binary link=logit solution;
	random int/ subject=ID;
	by sim;
run;

ods output ParameterEstimates=glimFH(keep = sim Effect Estimate) CovParms=glimCovH(keep = sim CovParm Estimate);
proc glimmix data=Sem method=LAPLACE;
	where ID<=500;
	class D_lag_0;
	model D_lag_0(event='1') = Y_lag_0 /dist=binary link=logit solution;
	random int/ subject=ID;
	by sim;
run;


DATA Glimcova;
set Glimcova;
scenario='A';
run;

DATA Glimcovb;
set Glimcovb;
scenario='B';
run;

DATA Glimcovc;
set Glimcovc;
scenario='C';
run;

DATA Glimcovd;
set Glimcovd;
scenario='D';
run;

DATA Glimcove;
set Glimcove;
scenario='E';
run;

DATA Glimcovf;
set Glimcovf;
scenario='F';
run;

DATA Glimcovg;
set Glimcovg;
scenario='G';
run;

DATA GlimcovH;
set GlimcovH;
scenario='H';
run;

DATA GlimcovI;
set Glimcov;
scenario='I';
run;

DATA Glimcov;
set Glimcova Glimcovb Glimcovc Glimcovd Glimcove Glimcovf Glimcovg Glimcovh Glimcovi;
run;

proc export data = Glimcov
outfile ="D:\Documents\PhD\Infering_idividual_causal\EX3_cov.csv"
dbms = csv
replace ;
run;

*fixed;
DATA Glimfa;
set Glimfa;
scenario='A';
run;

DATA Glimfb;
set Glimfb;
scenario='B';
run;

DATA Glimfc;
set Glimfc;
scenario='C';
run;

DATA Glimfd;
set Glimfd;
scenario='D';
run;

DATA Glimfe;
set Glimfe;
scenario='E';
run;

DATA Glimff;
set Glimff;
scenario='F';
run;

DATA Glimfg;
set Glimfg;
scenario='G';
run;

DATA GlimfH;
set GlimfH;
scenario='H';
run;

DATA GlimfI;
set Glimf;
scenario='I';
run;

DATA Glimf;
set Glimfa Glimfb Glimfc Glimfd Glimfe Glimff Glimfg Glimfh Glimfi;
run;

proc export data = Glimf
outfile ="D:\Documents\PhD\Infering_idividual_causal\EX3_fixed.csv"
dbms = csv
replace ;
run;

*Example 2;

%simdata(nsim=1, n=1000, nT=100, b0=0, b1=-0.2, b2=-0.1, bC=4, c0=0.1, c1=1, sigma=0.25, pC=0.5, tauC=0, tau0=0.25 , tau1=0.5 , tau2=0.25, tau3=0.25, tauCA1=0, tauCA2=0, a0=-0.5, a1=0.01, a2=1,aC=0.7); 

proc sort data=Sem;
by sim;
run;

data Ex2;
set SEM;
where sim=1;
run;

proc export data = Ex2
outfile ="D:\Documents\PhD\Infering_idividual_causal\SIM2.csv"
dbms = csv
replace ;
run;


ods output SolutionF=testF(keep = sim Effect Estimate) CovParms=testCov(keep = sim CovParm Estimate) SolutionR=testR(keep= sim Effect Subject Estimate  rename=(Effect=EffectR Estimate=EstimateR)) G=G V=V;
proc mixed data=SEM method=REML ASYCOV;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB outp=est_eps;;
	random int A_lag_2 A_lag_1 /subject=ID type=VC solution G V;
	by sim;
run;

ods output SolutionF=testFF(keep = sim Effect Estimate) CovParms=testCovF(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where ID<=500;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFC(keep = sim Effect Estimate) CovParms=testCovC(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where ID<=100;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFG(keep = sim Effect Estimate) CovParms=testCovG(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=3;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFD(keep = sim Effect Estimate) CovParms=testCovD(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=3 AND ID<=500;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFA(keep = sim Effect Estimate) CovParms=testCovA(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=3 AND ID<=100;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFH(keep = sim Effect Estimate) CovParms=testCovH(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=10;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFE(keep = sim Effect Estimate) CovParms=testCovE(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=10 AND ID<=500;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

ods output SolutionF=testFB(keep = sim Effect Estimate) CovParms=testCovB(keep = sim CovParm Estimate);
proc mixed data=SEM method=REML ASYCOV;
	where PER<=10 AND ID<=100;
	class PER;
	model Y_lag_0 = A_lag_2 A_lag_1 C_lag_1/solution COVB;
	random int A_lag_2 A_lag_1 /subject=ID type=VC;
	by sim;
run;

DATA Testcova;
set Testcova;
scenario='A';
run;

DATA Testcovb;
set Testcovb;
scenario='B';
run;

DATA Testcovc;
set Testcovc;
scenario='C';
run;

DATA Testcovd;
set Testcovd;
scenario='D';
run;

DATA Testcove;
set Testcove;
scenario='E';
run;

DATA Testcovf;
set Testcovf;
scenario='F';
run;

DATA Testcovg;
set Testcovg;
scenario='G';
run;

DATA TestcovH;
set TestcovH;
scenario='H';
run;

DATA TestcovI;
set Testcov;
scenario='I';
run;

DATA Testcov;
set Testcova Testcovb Testcovc Testcovd Testcove Testcovf Testcovg Testcovh Testcovi;
run;

proc export data = Testcov
outfile ="D:\Documents\PhD\Infering_idividual_causal\EX2_cov.csv"
dbms = csv
replace ;
run;

*fixed;
DATA Testfa;
set Testfa;
scenario='A';
run;

DATA Testfb;
set Testfb;
scenario='B';
run;

DATA Testfc;
set Testfc;
scenario='C';
run;

DATA Testfd;
set Testfd;
scenario='D';
run;

DATA Testfe;
set Testfe;
scenario='E';
run;

DATA Testff;
set Testff;
scenario='F';
run;

DATA Testfg;
set Testfg;
scenario='G';
run;

DATA TestfH;
set TestfH;
scenario='H';
run;

DATA TestfI;
set Testf;
scenario='I';
run;

DATA Testf;
set Testfa Testfb Testfc Testfd Testfe Testff Testfg Testfh Testfi;
run;

proc export data = Testf
outfile ="D:\Documents\PhD\Infering_idividual_causal\EX2_fixed.csv"
dbms = csv
replace ;
run;
