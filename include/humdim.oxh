#ifndef HDINCLUDED
	#define HDINCLUDED
	#include <oxstd.h> 
	#include <oxfloat.h>
	#import  <database>
	#import  <maximize>
	#include <oxprob.h>
	#include <oxdraw.h>
	#import  <solvenle>
	#import  <maxsqp>

//ENUMERATIONS

//1. Price vector enumeration
enum{equity,lab,res,permits,NPrices}

// Climate states
enum{atm1,atm2,atm3,Ntrans}
enum{stemp=Ntrans,otemp,OGHG,NClimate}	  //climate module

//Data set
enum{yeardat,gdpdat,emitdat,energydat,atmco2dat,NDatafiles}

//CF Dimensions of the problem
struct D {
	static decl
		A = 60,  // years of agent life
		P = 220, // periods
		Pm1= 219,
		Pm2= 218,
		G = 60   // generations
		;
		
	}
	
const decl CarbPrice = 50, C6K = 6000, C1m = 1E-3;

//Global functions

initialize();
parameterize();
exogenousfill(sc,pricefile,assetfile);
policy(permits_issued,carbon_tax,percap,sp,firm_allocation);
scenario(scene);
datacheck(j,NPP);
calibprint();
invest_nl(sysval,capital);
extract_nl(sysval,sent_vals);
agent_euler(sysval,startassets);
agents_problem(Ngen,itermax=1000);
newclimate(emissions);
equil(file_load="2",file_save="2_50");
parset(params);
moment_phi(params_sent,retval,score,hess);
moment(params_sent,retval,score,hess);
calibration();
caloutput();
polisim();

//Global Variables
const decl firm_toler = 1E-5,
			equil_toler = 1E-5,  //CF changed from -6, -7
			indir = "input\\",
			outdir = "output\\";	

struct Ext {
	static decl myX,xprice,kt,myout;
	static nl(sysval,sent_vals);
	}

struct Inv {
	static const decl kt_toler = 0.1;
	static decl kt, myK, mySH, myGDP, mpk, dividends;
	static nl(sysval,capital);
	static Production();
	}

struct GV {		//Generation-Specific Values
	static decl
		EE,
		LS,
		LS1,
		P,	//price matrix
		D,   // dividend stream
		A,	  // assetmatrix;
		PxA   // value matrix
		;
	static agent_euler(sysval,startassets);
	}
			

//Storage objects
decl data, params, prices;

//State Variables
decl SH,K,N,R,L,pop,popshare;
decl cstate;
decl CumC,mextcost;
decl lifespan; //CF ADDED this to avoid putting it in endow
decl assets,endow,hcap;
decl Omega,Omegat;
decl E,Mb,F,Ts;//environmental variables


//Policy variables
decl ctax,pstart,P,Pfirm;
decl rentinc,permitinc,dividends,shareprice,firm_permit_value;
decl rentshare,permitshare;


//Economic aggregates
decl U,GDP,cons,income,permit_transfers; //CF: took out utils
decl C,X;     //aggregate consumption and investment
decl resmarkup;//use this for calibration
decl Z; //number of firms in energy sector

//model parameters
decl rho,beta,sigma,delta,Ubar,sig1inv,del1inv; //utility function parameters
decl alpha,theta;
decl gammaa,gamman,rplus1;  //exogenous growth rates
decl deltaa,deltan;//growth decay rates
decl gammaphi,deltaphi;//growth decay rates
decl gammao,deltao;//growth decay rates
decl gammath,deltath,deltath2;//growth decay rates
decl phi;
decl t1,t2,t3,t4,deltam,t2x;
decl d0,b1,b2;
decl xi;


//solution objects

decl assetspass,pricespass,dividendspass,assetveclast;
decl newassets,euler;
decl limit_low,limit_high;
decl print_equil=0;
decl NP,NG;

#endif
