#import <humdim>

//ENUMERATIONS

//1. Price vector enumeration
enum{equity,lab,res,permits}
// Climate states
enum{atm1,atm2,atm3,stemp,otemp,OGHG}	  //climate module
//Data set
enum{yeardat,gdpdat,emitdat,energydat,atmco2dat}


//Global Variables


//Storage objects
extern decl data;
extern decl prices;  

//State Variables
extern decl K,N,R,L,pop,popshare;
extern decl cstate;
extern decl CumC,mextcost;
extern decl assets,endow,hcap;
extern decl Omega,Omegat;
extern decl E,Mb,F,Ts;//environmental variables
extern decl timet;  //to pass time period


//Policy variables
extern decl ctax,pstart,P,Pfirm;
extern decl rentinc,permitinc,dividends,shareprice,firm_permit_value;
extern decl rentshare,permitshare;


//Economic aggregates
extern decl U,GDP,cons,utils,income,permit_transfers;
extern decl C,X;     //aggregate consumption and investment
extern decl resmarkup;//use this for calibration
extern decl Z; //number of firms in energy sector

//model parameters
extern decl rho,beta,sigma,delta,Ubar; //utility function parameters
extern decl alpha,theta;
extern decl gammaa,gamman;  //exogenous growth rates
extern decl deltaa,deltan;//growth decay rates
extern decl gammaphi,deltaphi;//growth decay rates
extern decl gammao,deltao;//growth decay rates
extern decl gammath,deltath,deltath2;//growth decay rates
extern decl phi;
extern decl t1,t2,t3,t4,deltam,t2x;
extern decl d0,b1,b2;
extern decl xi;


//solution objects

extern decl assetspass,pricespass,dividendspass,assetveclast;
extern decl newassets,euler;
extern decl limit_low,limit_high;
extern decl print_equil;
extern decl NP,NG;

main(){
    decl starttime,fintime,params,pol,sc;
    decl ptest,j,over,i;
    starttime=timer();
    NG=60;  //number of generations
    NP=220*(NG/60);//number of periods
    pstart=48;
    decl dbase = new Database();
    dbase.LoadObs("edata.txt",5,60,0,1,1,1);
    data=dbase.GetAll();
    initialize(); //set up all variables
    parameterize();//parameterize the model
    exogenousfill(2,"./data/pricesfile.dat","./data/assetsfile.dat"); //fill vectors of exogenous state variables and policy parameters
    Z=2;//set level of market power in resource sector to acheive approx \$100 per ton markup

	calibration();
    println(deltam);
    fintime=timer();
    println("Total Run Time ",timespan(starttime,fintime));
    }
