#include "humdim.oxh"

//ENUMERATIONS

//1. Price vector enumeration
enum{equity,lab,res,permits}
// Climate states
enum{atm1,atm2,atm3,stemp,otemp,OGHG}	  //climate module
//Data set
enum{yeardat,gdpdat,emitdat,energydat,atmco2dat}


//Global Variables


//Storage objects
decl data;
decl prices;  

//State Variables
decl K,N,R,L,pop,popshare;
decl cstate;
decl CumC,mextcost;
decl assets,endow,hcap;
decl Omega,Omegat;
decl E,Mb,F,Ts;//environmental variables
decl timet;  //to pass time period


//Policy variables
decl ctax,pstart,P,Pfirm;
decl rentinc,permitinc,dividends,shareprice,firm_permit_value;
decl rentshare,permitshare;


//Economic aggregates
decl U,GDP,cons,utils,income,permit_transfers;
decl C,X;     //aggregate consumption and investment
decl resmarkup;//use this for calibration
decl Z; //number of firms in energy sector

//model parameters
decl rho,beta,sigma,delta,Ubar; //utility function parameters
decl alpha,theta;
decl gammaa,gamman;  //exogenous growth rates
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

//procedure to set up variables
initialize() {
  xi=new matrix[3];
  GDP=new matrix[NP];
  Omegat=new matrix[NP];
  U=new matrix[NP];
  prices=new matrix[NP][4];
  K=new matrix[NP];
  X=new matrix[NP];
  P=new matrix[NP];
  N=new matrix[NP][NG];
  pop=new matrix[NP][NG];
  cstate=new matrix[NP][6];
  deltam=new matrix[3][3];
  E=new matrix[NP];

  CumC=new matrix[NP];
  rho=new matrix[NP];
  assets=new matrix[NP][NG];
  cons=new matrix[NP][NG];
  utils=new matrix[NP][NG];
  newassets=new matrix[NP][NG];
  income=new matrix[NP]; //first period income storage
  permit_transfers=new matrix[NP];

  ctax=new matrix[NP];
  popshare=new matrix[NP][NG];//share of population represented by each agent in each cohort
  rentshare=new matrix[NP][NG];//share of resource rents
  permitshare=new matrix[NP][NG];//share of permit rents

  rentinc=new matrix[NP];//aggregate resource rents
  permitinc=new matrix[NP];//aggregate permit rents
  dividends=new matrix[NP];//total dividends (in equilibrium, number of shares sum to one, so it will be per share
  shareprice=new matrix[NP];//total dividends (in equilibrium, number of shares sum to one, so it will be per share

  phi=new matrix[NP];
  theta=new matrix[NP];
  }

  //procedure to assign parameters of the representative agent model
parameterize() {
decl i;
 beta=.96;
 delta=.045;//Pizer value
 Ubar=10;

 sigma=1.2213;	 //Pizer value
 alpha=.3;//Nordhaus Value

//United Nations Projections
 gamman=.02;
 deltan=.03;	//UN Medium Growth
 //deltan=0.018; //UN High Growth


    gammao=0.00762625; //calibrated
    deltao=1.07788e-006; //calibrated
    Omegat[]=0.0251492; //calibrated

    gammath=-0.0114053; //calibrated
    deltath=0.000550977; //calibrated
    theta[]=0.0553807; //calibrated

    deltath2=1;


    gammaphi=0.0117432; //calibrated
    deltaphi=0.0509572; //calibrated
    phi[0]=1.30726; //calibrated


    t2x=2.98;
    t1=.9112;//Pizer Value
    t3=0.002; //From Dice99, convert to ten year
    t4=.01012;	 //from Dice99
    t2x=2.98;
    t2=t2x*(1-t4-t1)/log(2);//pizer and nordhaus
    //Damage Coefficients
    b1=-0.004500;
    b2=0.003500;
    d0=1;

    deltam[0][0]=exp(log(.66616)/10);
    deltam[0][1]=1-deltam[0][0];
    deltam[1][1]=exp(log(.60897)/10);
    deltam[1][0]=.27607/(.27607+.11496)*(1-deltam[1][1]);
    deltam[1][2]=.11496/(.27607+.11496)*(1-deltam[1][1]);
    deltam[2][2]=exp(log(.99578)/10);
    deltam[2][1]=1-deltam[2][2];

    xi[0]=113;
    xi[1]=700;
    xi[2]=4;

    //load earnings profile from file profile.dta
    decl db,hcbase;
    db = new Database();
    if(db.LoadDta("./data/profile.dta",1,1,1)==FALSE){
		oxrunerror("Could not read: ./data/profile.dta", 0);	
	}
    hcbase=db.GetAll();
    hcap=hcbase[][1];
	
    }

 //given values, fill in vectors for At, L
exogenousfill(sc,pricefile,assetfile) {
  decl i,ii;

 //Initial Values
    CumC[]=213;
    Mb=590;
    cstate[][atm1]= 690;
    cstate[][atm2]= 770;
    cstate[][atm3]= 19180;
    cstate[][otemp]=0;
    cstate[][stemp]=0.295;
    rho[]=1;

  for(i=1;i<NP;++i) {
    phi[i]=phi[i-1].*(1+gammaphi*(1-deltaphi)^(i));
	Omegat[i]=Omegat[i-1].*(1+gammao*(1-deltao)^(i));
	theta[i]=theta[i-1]*(1+gammath*(1-deltath*deltath2^i)^i);
	}
Omega=Omegat;//set initial Omega to trend value - bring in damages in equil()
decl j;
pop[0][0]=88.5*10^6;
for(j=1;j<NG;++j)
	pop[0][j]=pop[0][j-1]./(1+gamman);
for(i=1;i<NP;++i)
    {
	pop[i][0]=pop[i-1][0].*(1+gamman*(1-deltan)^(i));
        for(j=1;j<NG;++j)
        pop[i][j]=pop[i-1][j-1];
    }
N=pop.*hcap';//N is units of effective labour supply per cohort per time period
L=sumr(N);  //L is effective labour supply in a time period
N=N./10^6; //labour supply by cohort
L=L./10^6; //aggregate labour supply.
pop=pop./10^6;

popshare=ones(NP,NG)./sumr(pop.*10^6);  //each cohort's share of the total population
//each agent in a particular cohort makes up popshare of the population
//cohort share of pop is popshare.*pop.*10^6
rentshare=popshare;
permitshare=popshare;
prices[][equity]=constant(1,NP,1);//start with a simple interest-free savings
//prices[15:NP-1][equity]=prices[15][equity];
prices[][lab]=.004;
prices[][res]=113.5;
prices[][permits]=0;

	decl dbase;
	dbase = new Database();
	if(dbase.Load(assetfile)==FALSE){
		oxrunerror("Could not read: "~assetfile, 0);	
	}
	assets=dbase.GetAll();
	dbase = new Database();
	if(dbase.Load(pricefile)==FALSE){
		oxrunerror("Could not read: "~pricefile, 0);	
	}
	prices=dbase.GetAll();
	return 1;
}

policy(permits_issued,carbon_tax,percap,sp,firm_allocation) {
	println("in policy ",permits_issued~carbon_tax~percap~sp~firm_allocation);

	P[]=50;//set initial permits to be non-binding
	P[sp:NP-1]=permits_issued;
	if(permits_issued==50) {
		//println("load 50 dollar shadow price file");
		decl dbase = new Database();
		if(dbase.Load("./data/shadow_50.dat")){
			oxrunerror("Could not read: ./data/shadow_50.dat", 0);	
		}
		P[sp:NP-1]=dbase.GetAll()[sp:NP-1];
		}	

	if(percap==1)
		P[sp:NP-1]=P[sp:NP-1]./sumr(pop[sp]).*sumr(pop[sp:NP-1]);  //per capita permits

	ctax[]=0;//the carbon tax rate is zero for early periods
	ctax[sp:NP-1]=carbon_tax;//set carbon tax for policy period
	Pfirm=P.*firm_allocation;

	}

scenario(scene) {
	switch_single(scene) {
	case 1: //low sensitivity
		t2x=1.5;
    	t2=t2x*(1-t4-t1)/log(2);
    	d0=1;
	case 2:  //pizer and nordhaus
		t2x=2.98;
    	t2=t2x*(1-t4-t1)/log(2);
    	d0=1;
	case 3: 	//high sensitivity
    	t2x=4.5;
		t2=t2x*(1-t4-t1)/log(2);
    	d0=1;
    }
}

datacheck(j,NPP) {
	decl i,errors;
	errors=new matrix[60][3];
	for(j=0;j<60;++j)	{
		i=j+9;
		println(1962+i,"  ",GDP[i],"    ",data[j][gdpdat],"    ",E[i]*phi[i],"   ",data[j][energydat],"   ",E[i],"   ",data[j][emitdat],"   ",cstate[i][atm1],"   ",data[j][atmco2dat]);
		errors[j][0]=(GDP[i]-data[j][gdpdat])./data[0][gdpdat];
		errors[j][1]=(phi[i]*E[i]-data[j][energydat])./data[0][energydat];//get to same order of magnitude as GWP - use first period levels
		errors[j][2]=(phi[i]-data[j][energydat]/data[j][emitdat])./(data[0][energydat]/data[0][emitdat]);
		//errors[j][1:2].*=1;
		}
	return errors[][];//return scaled errors for energy and gwp
	}


calibprint() {
	decl i;
	for(i=1;i<NP;++i)
		println(i,"  ",Omegat[i],"   ",Omega[i],"   ",cstate[i][],"     ",E[i],"     ",CumC[i],"     ",mextcost[i],"     ",E[i]*phi[i],"       ",L[i],"      ",double(sumr(pop[i][])),"      ",pop[i][0],"   ",K[i],"   ",U[i],"    ",GDP[i],"    ",cons[i][0],"   ",sumr(cons[i][].*pop[i][0]*10^(-6))[0][0],"   ",phi[i][0],"   ",theta[i][0],"   ",prices[i][0],"   ",prices[i][1],"   ",prices[i][2],"   ",prices[i][3],"   ",P[i][0],"   ",Pfirm[i][0],"   ",ctax[i][0],"   ",rentinc[i][0],"   ",permitinc[i],"   ",dividends[i],"   ",income[i],"   ",permit_transfers[i]);
		}

 //evaluate system of agent's euler equations
invest_nl(sysval,capital) {
	decl euler,surplus,kt;
	decl sent=shape(capital,NP-1,2);
	decl capchoice=(K[0]|sent[][0])~(0|sent[][1]);
	decl mpk=alpha.*Omega.*capchoice[][0].^(alpha-1).*L.^(1-alpha-theta).*(phi.*X).^(theta)+ones(NP,1)-capchoice[][1];
	decl discount=beta.*ones(NP,1)./prices[][equity];

	GDP=Omega.*capchoice[][0].^(alpha).*L.^(1-alpha-theta).*(phi.*X).^(theta);
	prices[][lab]=(1-alpha-theta).*GDP./L;
	prices[][res]=(theta).*GDP./X;

	surplus=Omega.*capchoice[][0].^(alpha).*L.^(1-alpha-theta).*(phi.*X).^(theta)-prices[][lab].*L-prices[][res].*X+prices[][permits].*Pfirm+capchoice[][0];
	dividends=surplus-lag0(capchoice[][0],-1)/(1-delta);//dividends per share from production in dollars
	dividends[NP-1] -= capchoice[][0][NP-1][0]./(1-delta); //force firm to maintain capital stock.

	//use .1 here - as long as it is not binding in equil it is fine...using exactly zero admits small neg values.
	kt=(dividends[:NP-2].<=0.1).*(dividends[:NP-2]-0.1)+(dividends[:NP-2].>0.1).*capchoice[1:][1]+(capchoice[1:][1].>0).*capchoice[1:][1];
	euler=discount[:NP-2].*mpk[1:][]-ones(NP-1,1)./(1-delta);

	if(print_equil==1)
		println((0|euler)~(dividends)~capchoice~mpk~X);
	sysval[0]=(euler|kt);
	}

 //evaluate system of agent's euler equations
extract_nl(sysval,sent_vals) {
	decl i,mpr,revenue;
	decl sent=shape(sent_vals,NP,2);
	decl extraction=exp(sent[][0]);//constrain to be positive
	decl Cumtemp=cumsum(lag0(extraction,1),1);
	decl Cum_eps=cumsum(lag0(extraction+((10^(-3))|zeros(NP-1,1)),1),1);

	decl costs=10^(-3)*xi[0]*extraction+10^(-3)*(xi[1]*6000/(xi[2]+1))*(((Cumtemp+extraction)/6000).^(xi[2]+1)-((Cumtemp)/6000).^(xi[2]+1));
	//calculate the change in extraction costs for a small increase in extraction
	decl costs_eps=10^(-3)*xi[0]*extraction+10^(-3)*(xi[1]*6000/(xi[2]+1))*(((Cum_eps+extraction)/6000).^(xi[2]+1)-((Cum_eps)/6000).^(xi[2]+1));

	decl costinc=reversec(cumsum(reversec(costs_eps-costs),ones(NP,1)*beta./reversec(prices[][equity])))./10^(-3);

	decl opp_cost=costinc;
	//opp_cost is the current value of all future cost increases caused by a one unit change in extraction in that period

	decl marg_rev=(theta.*Omega.*K.^alpha.*L.^(1-alpha-theta).*(phi.*extraction).^(theta)./extraction)+1/Z*extraction.*(theta.*(theta-1).*Omega.*K.^alpha.*L.^(1-alpha-theta).*(phi.*extraction).^(theta)./(extraction.^2));//mkt power - prices are endogenous
	marg_rev+=sent[][1]-ctax;
	mextcost=10^-3*(xi[0]+xi[1].*((Cumtemp[][0]+extraction)/6000).^xi[2]);//marginal extraction cost

	decl euler=marg_rev-mextcost-1/Z*opp_cost;

	decl kt=(extraction.>P).*(extraction-P).^2+(extraction.<P).*sent[][1].^2+(sent[][1].>0).*sent[][1].^2;
	prices[][permits]=-sent[][1];

	mpr=theta.*Omega.*K.^(alpha).*L.^(1-alpha-theta).*(phi.*X).^(theta)./X;
	prices[][res]=mpr-prices[][permits]-ctax;
	resmarkup=mpr-prices[][permits]-ctax-mextcost;
	revenue=prices[][res].*extraction;
	rentinc=(revenue-costs).*10^12;//total rent income in dollars to pass to agents' problem
	sysval[0]=euler|kt;
	if(print_equil==1)
		println("resource sector",extraction~prices[][res]~P~prices[][permits]~rentinc);
	}

 //evaluate system of agent's euler equations
new_assetsnl(sysval,startassets) {
	//enum{equity,lab,res,permits}
	decl lifespan=endow[0][0];
	decl assetstemp=((endow[0][1]|startassets)~assetspass[][1:3])|zeros(1,4);//first period endowment and death period zero added
	decl pricestemp=pricespass|pricespass[rows(pricespass)-1][];
	decl divstemp=(dividendspass|0).*assetstemp[][0];//dividends per share times shareholdings

	decl consumption=sumr(pricestemp.*assetstemp)-pricestemp[][equity].*lag0(assetstemp[][0],-1)+divstemp[][0];
	sysval[0]=(consumption[0:lifespan-2][]).^(-sigma)-beta./pricestemp[0:lifespan-2][equity].*(pricestemp[1:lifespan-1][equity]+dividendspass[1:lifespan-1][]).*(consumption[1:lifespan-1][]).^(-sigma);
	sysval[0]*=10;//scale it for solution tightness
	//println("relevant measures - consumption, income, investment, dividends, euler",consumption~(sumr(pricestemp.*assetstemp))~(pricestemp[][equity].*lag0(assetstemp[][0],-1))~divstemp[][0]~(0|sysval[0]|0));
}
						
//equity distribution iteration
agents_problem(itermax)	{
	decl fstart;
	decl born,lifespan,startlife,converge,euler;
	decl storevec,assetvec,assetfill,consumption,consfill,utilityfill;
	decl test,np_age,sp_ext,step;
	converge=0;
	//println("Beginning agent's problem");
	sp_ext=new matrix[NG][1];//extra NG periods of share prices, holding last period ROR and dividend constant
	sp_ext[0][0]=shareprice[NP-1]*prices[NP-1][equity]-dividends[NP-1];
	for(step=1;step<NG;++step)  {
	  sp_ext[step][0]=sp_ext[step-1]*prices[NP-1][equity]-dividends[NP-1];
 	  }//end for loop													

	//println(shareprice|sp_ext);

	for(born=-NG+1;born<NP;++born) {
		fstart=timer();
		lifespan=NG;
		startlife=0;
		np_age=60;
		if(born<0) {
			startlife=-born;
			lifespan=NG-startlife;
			}
		if(born<=NP-NG)	{
			//println("  ");
			assetspass=(diagonal(assets[max(born,0):NP-1][startlife:NG-1])')~hcap[startlife:startlife+lifespan-1][0]~diagonal(rentshare[max(born,0):NP-1][startlife:NG-1])'~diagonal(permitshare[max(born,0):NP-1][startlife:NG-1])';
			//prices paid on assets are equity price+per share dividend, wage rate, total res rents, total permit rents
			pricespass=(shareprice[max(born,0):max(born,0)+lifespan-1][0])~(prices[max(born,0):max(born,0)+lifespan-1][lab]*10^6)~(rentinc[max(born,0):max(born,0)+lifespan-1][0])~permitinc[max(born,0):max(born,0)+lifespan-1][0];			
			dividendspass=dividends[max(born,0):max(born,0)+lifespan-1][0];//pass on relevant time period dividends per share
			}
		if(born==NP-NG-2) storevec=assetspass;
		if(born>NP-NG){
			startlife=0;
			lifespan=NG;
			np_age=NP-born;
			assetspass=((diagonal(assets[max(born,0):NP-1][startlife:NG-1])')|zeros(NG-np_age,1))~hcap[startlife:startlife+lifespan-1][0]~((diagonal(rentshare[max(born,0):NP-1][startlife:np_age-1])')|(ones(NG-np_age,1).*rentshare[NP-1][np_age-1]))~((diagonal(permitshare[max(born,0):NP-1][startlife:np_age-1])')|(ones(NG-np_age,1).*permitshare[NP-1][np_age-1]));
			assetspass[][0]=storevec./1000;//try these as starting values - should be okay.
			pricespass=(shareprice[max(born,0):min(max(born,0)+lifespan-1,NP-1)][0]|(sp_ext[0:NG-np_age-1][0]))~((prices[max(born,0):min(max(born,0)+lifespan-1,NP-1)][lab]*10^6)|(prices[NP-1][lab]*10^6 .*ones(NG-np_age,1)))~(rentinc[max(born,0):min(max(born,0)+lifespan-1,NP-1)][0]|(rentinc[NP-1][0].*ones(NG-np_age,1)))~(permitinc[max(born,0):min(max(born,0)+lifespan-1,NP-1)][0]|(permitinc[NP-1][0].*ones(NG-np_age,1)));			
			dividendspass=dividends[max(born,0):min(max(born,0)+lifespan-1,NP-1)][0]|(dividends[NP-1][0].*ones(NG-np_age,1));//pass on relevant time period dividends per share with fixed last period dividend through time
			}				
		//println("Born in period ",born," with lifespan ",lifespan);
		//println("age at NP=",np_age," padded with ",NG-np_age,"zeros.");
		//assets are equity, human capital, res rent allocation, and permit rent allocation
		consumption=sumr(pricespass.*assetspass);
		assetvec=vec(assetspass[0][0]);
		if(lifespan==1) {
			cons[max(born,0)][startlife]=consumption[0][0];
			utils[max(born,0)][startlife]=consumption[0][0].^(1-sigma)./(1-sigma)+Ubar;
			//println("Consumption ",(cumsum(ones(rows(consumption),1),1)-1)~(sumr(pricespass.*assetspass))~consumption);
			}
		else  {  //CF: changed to else from "if(lifespan>1)"
			endow=lifespan~assetspass[0][];//row vector of all three endowments
			assetvec=vec(assetspass[0:lifespan-1][0]);//vector of capital holdings, in partial shares
			test=assetvec;
			MaxControl(1000,0);
			euler=assetvec;
			converge+=(SolveNLE(new_assetsnl,&assetvec));
			new_assetsnl(&euler,assetvec);
			assetspass[][0]=((endow[0][1]|assetvec));
			consumption=sumr(pricespass.*assetspass)-pricespass[][equity].*lag0(assetspass[][0],-1)+(dividendspass.*assetspass[][0]);
			assetfill=assets[max(born,0):min(max(born,0)+lifespan-1,NP-1)][startlife:];
			consfill=cons[max(born,0):min(max(born,0)+lifespan-1,NP-1)][startlife:];
			assetfill=setdiagonal(assetfill,endow[0][1]|submat(assetvec,0,np_age-2,0,0));
			consfill=setdiagonal(consfill,submat(consumption,0,np_age-1,0,0));

			assets[max(born,0):min(max(born,0)+lifespan-1,NP-1)][startlife:]=assetfill;
			cons[max(born,0):min(max(born,0)+lifespan-1,NP-1)][startlife:]=consfill;

		if(born>=0) {
			U[born]=sumc((consumption.^(1-sigma)./(1-sigma)+Ubar).*.95.^(cumsum(ones(rows(consumption),1),1)-1));
			permit_transfers[born][0]=sumr(pricespass[0][3].*assetspass[0][3]);
			income[born][0]=sumr(pricespass[0][].*assetspass[0][])+(dividendspass[0][].*assetspass[0][0]);
			}
		//println("  Convergence flag from agent's problem:",converge);
		}
	}//end of for loop over birth cohorts

	utils=cons.^(1-sigma)./(1-sigma)+Ubar;

	println("Sum of Convergence flags from agent's problem (must be zero):",converge);
	return converge;
	}

// Evolution of the climate model given emissions.
newclimate(emissions)	{
	//use global NPx6 matrix of climate values and receive a NPx1 vector of emissions levels
	//fill a new NPx6 matrix of updated climate values.
	decl i;
	//cstate-=<596.4,705,19200,0,0,0>;//convert to deviations from P.I. normals
	for(i=1;i<NP;++i)	{
		cstate[i][atm1]=emissions[i-1][0]+deltam[atm1][atm1]*cstate[i-1][atm1]+deltam[atm2][atm1]*cstate[i-1][atm2];
		cstate[i][atm2]=deltam[atm2][atm2]*cstate[i-1][atm2]+deltam[atm1][atm2]*cstate[i-1][atm1]+deltam[atm3][atm2]*cstate[i-1][atm3];
		cstate[i][atm3]=deltam[atm3][atm3]*cstate[i-1][atm3]+deltam[atm2][atm3]*cstate[i-1][atm2];
		//cstate+=<596.4,705,19200,0,0,0>;  //transform back to levels
		//c_new=c_new;//
		cstate[i][stemp]=t4*cstate[i-1][otemp]+t1*cstate[i-1][stemp]+t2*log((cstate[i-1][atm1])/596.4);
		cstate[i][otemp]=cstate[i-1][otemp]+t3*(cstate[i-1][stemp]-cstate[i-1][otemp]);
		}
	}


//main equilibrium procedure
equil(file_load,file_save) {
	decl crit_equil,agent_crit,firm_crit,agent_converge;
	decl kold,xold,adj,ror_old;
	decl extstart,extchoice,capstart,capchoice;
	decl extcost,extrevenue;
	decl mpk,mpl,mpr;
	decl dbase;
	decl mcost;
	decl K_converge=0;
	decl X_converge=0;

	//begin loop over asset rate of return/shareholder discount factor.
	K=19/1.005*cumprod(constant(1.005,NP,1));//initial guess of capital stock increases at the rate of discount
	X=2*cumprod(constant(1.005,NP,1));//initial guess of capital stock increases at the rate of discount
	decl ext_euler=X;
	decl inv_euler=K;

	dbase = new Database();
	if(dbase.Load("./data/Kfile"~file_load~".dat")==FALSE){
		oxrunerror("Could not read: ./data/Kfile"~file_load~".dat", 0);	
	}
	K=dbase.GetAll();

	dbase = new Database();
	if(dbase.Load("./data/Xfile"~file_load~".dat")==FALSE){
		oxrunerror("Could not read: ./data/Xfile"~file_load~".dat", 0);	
	}
	X=dbase.GetAll();

	dbase = new Database();
	if(dbase.Load("./data/assetsfile"~file_load~".dat")==FALSE){
		oxrunerror("Could not read: ./data/assetsfile"~file_load~".dat", 0);	
	}
	assets=dbase.GetAll();

	dbase = new Database();
	if(dbase.Load("./data/pricesfile"~file_load~".dat")==FALSE){
		oxrunerror("Could not read: ./data/pricesfile"~file_load~".dat", 0);	
	}
	prices=dbase.GetAll();

	Omega=Omegat;

	crit_equil=10;
	while(crit_equil>10^(-7))	{
		MaxControl(200,0);
		ror_old=prices[][equity];
		firm_crit=7;
		//solve extraction and production given rate of return
		while(firm_crit>10^(-6)) {									
			kold=K;
			xold=X;
			E=X;
			newclimate(E);//send sequence of emissions to temp routine
			//update climate damages
			Omega=Omegat./(constant(d0,NP,1).*(ones(NP,1)+b1*cstate[][stemp]+b2*cstate[][stemp].^2));
			extstart=log(X)|zeros(NP,1);
			//println("Calculated Res Ext, given K with convergence crit ", SolveNLE(extract_nl,&extstart,-1));
        	X_converge=SolveNLE(extract_nl,&extstart,-1);
			print_equil=0;
			extract_nl(&ext_euler,extstart);
			print_equil=0;
			extchoice=shape(extstart,NP,2);
			X=exp(extchoice[][0]);
			CumC=cumsum(lag0(X,1),1);
			extrevenue=(theta.*Omega.*K.^alpha.*L.^(1-alpha-theta).*(phi.*X).^(theta));
			extcost=10^(-3)*xi[0]*X+10^(-3)*(xi[1]*6000/(xi[2]+1))*(((CumC+X)/6000).^(xi[2]+1)-((CumC)/6000).^(xi[2]+1));

			capstart=(K[1:NP-1][])|zeros(NP-1,1);//capital and shadow values
			//println("Calculated K, given res ext with convergence crit ",K_converge=SolveNLE(invest_nl,&capstart,-1));
			K_converge=SolveNLE(invest_nl,&capstart,-1);
			print_equil=0;
			invest_nl(&inv_euler,capstart);
			print_equil=0;
			capchoice=shape(capstart,NP-1,2);
			K=(K[0][0]|(capchoice[][0]));
		
			//update omega
			//println("Updating climate given res ext");
			prices[][permits]=-extchoice[][1];
			firm_crit=meanc(meanr(((K~X)-(kold~xold))).^2); //mean square deviation
			xold=X;
			kold=K;
			println("Firm Convergence crit, X, K ",firm_crit[0][0],"   ",X_converge[0][0],"   ",K_converge[0][0]);
			}
		print_equil=0;
		extract_nl(&ext_euler,extstart);
		print_equil=0;
		permitinc=(P!=Pfirm).*(X-Pfirm).*prices[][permits]*10^12;
		firm_permit_value=-(Pfirm).*prices[][permits];
		//update share prices given dividends and rates of return
		shareprice[][]=1000;//the price of one share in the last period in dollars
		for(timet=NP-2;timet>=0;--timet)  {
		  shareprice[timet][0]=(shareprice[timet+1][0]+dividends[timet+1])/prices[timet+1][equity];
	 	  }//end for loop
		agent_converge=0;
		agent_crit=10^(-8);
		agent_converge=agents_problem(1000);
		while(agent_converge>0) {
			println("Unconverged assets - recalculating");
			assets[1:][]=0;
     		agent_converge=agents_problem(1000);
			}
		adj=((sumr(assets.*pop.*10^6))./(10^12))-ones(NP,1);
		adj[][]=spline(adj,cumsum(ones(NP,1),1),0);//locally smooth adjustment - keeps it from going pos/neg/pos
		agent_crit=maxc(fabs(adj[][]));
		prices[][equity]./=(ones(NP,1)+adj[][]/50);//adjust rates of return
		if(savemat("./data/xfile"~file_save~".dat",X)==0){
			oxrunerror("Could not write to: ./data/xfile.dat", 0);	
		}
		if(savemat("./data/kfile"~file_save~".dat",K)==0){
			oxrunerror("Could not write to: ./data/kfile.dat", 0);	
		}
		if(savemat("./data/assetsfile"~file_save~".dat",assets)==0){
			oxrunerror("Could not write to: ./data/xfile.dat", 0);	
		}
		if(savemat("./data/pricesfile"~file_save~".dat",prices)==0){
			oxrunerror("Could not write to: ./data/xfile.dat", 0);	
		}
		crit_equil=meanc(fabs((prices[][equity])-(ror_old))[:NP-3][])+agent_crit; //firm crit will be small due to while loop above
		println("equil_crit and conv flags",crit_equil,"  X=",X_converge,"   K=",K_converge,"  agent=",agent_converge);
		}
	}



//procedure to assign parameters of the model for empirical calibration
parset(params) {
	gammao=params[0][0];
	deltao=params[1][0];
	Omegat[][]=params[2][0];

	gammath=params[3][0];
	deltath=params[4][0];
	theta[][]=params[5][0];

	gammaphi=params[6][0];
	deltaphi=params[7][0];
	phi[0][0]=params[8][0];

	exogenousfill(2,"./data/pricesfile.dat","./data/assetsfile.dat"); //fill vectors of exogenous state variables and policy parameters

	}

moment_phi(params_sent,retval,score,hess) {
	decl test;
	decl params=gammao|deltao|Omegat[0][0]|gammath|deltath|theta[0][0]|params_sent;
	parset(params);

	//equil();	//see if this works...you may need to fill E, GDP, etc to make it work.
	test=datacheck(0,0);

	fopen("calibration.log","l");
	println("Low limit   Estimates   High Limit");
	println(limit_low~params~limit_high);
	retval[0]=-sumc(sumr(test[][2].^2));
	println("Mean absolute value percentage deviation ",retval[0]);
	println(" //Calibrated Parameters - Output");
	println(" gammao=",params[0][0],"; //calibrated");
	println(" deltao=",params[1][0],"; //calibrated");
	println("Omegat[][]=",params[2][0],"; //calibrated");

	println(" ");
	println("gammaphi=",params[6][0],"; //calibrated");
	println("deltaphi=",params[7][0],"; //calibrated");
	println(" phi[0][0]=",params[8][0],"; //calibrated");
	println(" ");
	println(" ");
	println("gammath=",params[3][0],"; //calibrated");
	println("deltath=",params[4][0],"; //calibrated");
	println("theta[][]=",params[5][0],"; //calibrated");
	println(" ");
	println("deltath2=1;");
	fclose("l");
	//limits used to keep parameters in possible ranges - positive TFP, theta, and emissions ratios, positive or negative initial growth rates, etc.

	}

moment(params_sent,retval,score,hess) {
	decl test;
	decl params=params_sent|gammaphi|deltaphi|phi[0][0];
	parset(params);
	equil("2","calib");
	test=datacheck(0,0);
	fopen("calibration.log","l");
	println("Low limit   Estimates   High Limit");
	println(limit_low~params~limit_high);
	retval[0]=-sumc(sumr(test[][0:1].^2));
	println("Mean absolute value percentage deviation ",retval[0]);
	println(" //Calibrated Parameters - Output");
	println(" gammao=",params[0][0],"; //calibrated");
	println(" deltao=",params[1][0],"; //calibrated");
	println("Omegat[][]=",params[2][0],"; //calibrated");

	println(" ");
	println("gammaphi=",params[6][0],"; //calibrated");
	println("deltaphi=",params[7][0],"; //calibrated");
	println(" phi[0][0]=",params[8][0],"; //calibrated");
	println(" ");
	println(" ");
	println("gammath=",params[3][0],"; //calibrated");
	println("deltath=",params[4][0],"; //calibrated");
	println("theta[][]=",params[5][0],"; //calibrated");
	println(" ");
	println("deltath2=1;");
	fclose("l");
	//limits used to keep parameters in possible ranges - positive TFP, theta, and emissions ratios, positive or negative initial growth rates, etc.
	}

calibration() {
	decl params=gammao|deltao|Omega[0][0]|gammath|deltath|theta[0][0];
	decl params_phi=gammaphi|deltaphi|phi[0][0];
	limit_low=(params|params_phi).*.9;
	limit_low[3][0]=gammath*1.1;//gammath is negative
	limit_high=(params|params_phi).*1.1;
	limit_high[3][0]=gammath*.9;//gammath is negative

	decl err_val;
	policy(54,0,0,48,0);
	scenario(2);

	//MaxSQP(moment_phi,&params_phi,&err_val, 0, 1,0,0, limit_low[6:8][],limit_high[6:8][]);


	MaxSQP(moment,&params,&err_val, 0, 1,0,0, limit_low[0:5][],limit_high[0:5][]);
	//
	////savemat("paramsfile.dat",params);

	params=params|params_phi;

	println(" //Calibrated Parameters - Output");
	println(" ");
	println(" gammao=",params[0][0],"; //calibrated");
	println(" deltao=",params[1][0],"; //calibrated");
	println("Omegat[][]=",params[2][0],"; //calibrated");
	println(" ");
	println(" ");
	println("gammath=",params[3][0],"; //calibrated");
	println("deltath=",params[4][0],"; //calibrated");
	println("theta[][]=",params[5][0],"; //calibrated");
	println(" ");
	println("deltath2=1;");
	println(" ");
	println(" ");
	println("gammaphi=",params[6][0],"; //calibrated");
	println("deltaphi=",params[7][0],"; //calibrated");
	println(" phi[0][0]=",params[8][0],"; //calibrated");
	println(" ");
	println(" ");

}

//  Calibrated model output
caloutput() {


	policy(54,0,0,48,0);
	scenario(2);
	equil("","");
	fopen("calib_rev2.log","l");
	calibprint();
	fclose("l");
	if(savemat("./data/xfile.dat",X)==0){
		oxrunerror("Could not write to: ./data/xfile.dat", 0);	
	}
	if(savemat("./data/kfile.dat",K)==0){
		oxrunerror("Could not write to: ./data/kfile.dat", 0);	
	}
	if(savemat("./data/assetsfile.dat",assets)==0){
		oxrunerror("Could not write to: ./data/assetsfile.dat", 0);	
	}
	if(savemat("./data/pricesfile.dat",prices)==0){
		oxrunerror("Could not write to: ./data/pricesfile.dat", 0);	
	}

	}

//Policy Simulations
polisim() {	

	decl logfile,g,h,i,j;

	////policy(permits_issued,carbon_tax,percap,sp,firm_allocation)
	decl quotachoices=<540;54;500>;
	decl percapchoices=<0;1>;
	decl firmchoices=<0;1>;
	decl scchoices=<2>;


	//find \$50 shadow value emissions
	policy(54,0.0485,0,48,0);
	scenario(2);
	equil("2","2_50");
	if(savemat("./data/shadow_50.dat",X)==0){
		oxrunerror("Could not write to: ./data/xfile.dat", 0);	
	}

	println(ctax~prices[][permits]~X);

	policy(5.4,0,0,48,0);
	scenario(2);
	equil("","");
	if(savemat("./data/shadow_50.dat",X)==0){
		oxrunerror("Could not write to: ./data/xfile.dat", 0);	
	}

	println(ctax~prices[][permits]~X);

	for(g=2;g<rows(quotachoices);++g)
		for(j=0;j<rows(scchoices);++j)
			{
			if(g>0)
				{
				for(h=0;h<rows(percapchoices);++h)
					for(i=0;i<rows(firmchoices);++i)
						{
						//println("permit allocation ",quotachoices[g]," per capita indication ",percapchoices[h][0]," firm permits ",firmchoices[i][0], " and climate scenario ",scchoices[j][0]);
						policy(quotachoices[g][0]/10,0,percapchoices[h][0],48,firmchoices[i][0]);
						//println("permit allocation and firm permits",P~Pfirm);
						logfile="simuls_"~sprint(quotachoices[g][0])~"_"~sprint(percapchoices[h][0])~"_"~sprint(firmchoices[i][0])~"_"~sprint(scchoices[j][0])~".log";
						println("logfile is ",logfile);
						scenario(scchoices[j][0]);
						println("t2x and temp change parameter ",t2x~t2);
						equil(sprint(scchoices[j][0]),sprint(scchoices[j][0])~"_"~sprint(quotachoices[g][0]));
						fopen(logfile,"l");
						calibprint();
						fclose("l");
						}
				}
			else
				{
//					//println("permit allocation ",quotachoices[g]," per capita indication ",percapchoices[h][0]," firm permits ",firmchoices[i][0], " and climate scenario ",scchoices[j][0]);
				policy(54,0,0,48,0);
//					//println("permit allocation and firm permits",P~Pfirm);
				logfile="simuls_"~sprint(quotachoices[g][0])~"_"~sprint(0)~"_"~sprint(0)~"_"~sprint(scchoices[j][0])~".log";
				println("logfile is ",logfile);
				scenario(scchoices[j][0]);
				println("t2x and temp change parameter ",t2x~t2);
				equil(sprint(scchoices[j][0]),sprint(scchoices[j][0]));
				fopen(logfile,"l");
				calibprint();
				fclose("l");
				}
			}
		

		fopen("cohorts.log","l");
		println(cumsum(ones(NP,1),1)~pop[][0]);
		fclose("l");
	}
