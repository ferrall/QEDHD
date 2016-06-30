#include "humdim.oxh"

//procedure to set up variables
initialize() {
  xi=new matrix[3];

  prices=new matrix[D::P][4];

  E=  phi=  theta=
 income= permit_transfers= ctax=CumC= rho= GDP = Omegat=U= K= X= P= zeros(D::P,1);

  
  deltam=new matrix[Ntrans][Ntrans];
						   //CF: not needed?? utils=  
  N= pop=  assets=  cons=  newassets=
  popshare=  //share of population represented by each agent in each cohort
  rentshare= //share of resource rents
  permitshare= //share of permit rents
  	zeros(D::P,D::G);

  }

  //procedure to assign parameters of the representative agent model
parameterize() {
    decl i;
    beta=.96;
	rplus1 = 1.005;
    delta=.045;//Pizer value
	del1inv = 1/(1-delta); //CF added
	
    Ubar=10;

    sigma=1.2213;	 //Pizer value
	sig1inv = 1/(1-sigma);  //CF added 
    alpha=.3;//Nordhaus Value

    //United Nations Projections
    gamman=.02;
    deltan=.03;	//UN Medium Growth
    //deltan=0.018; //UN High Growth

    gammao=0.00762625; //calibrated
    deltao=1.07788E-006; //calibrated
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
	
    deltam[atm1][atm1]=exp(log(.66616)/10);
    deltam[atm2][atm1]=1-deltam[atm1][atm1];
    deltam[atm2][atm2]=exp(log(.60897)/10);
    deltam[atm2][atm1]=.27607/(.27607+.11496)*(1-deltam[atm2][atm2]);
    deltam[atm2][atm3]=.11496/(.27607+.11496)*(1-deltam[atm2][atm2]);	
    deltam[atm3][atm3]=exp(log(.99578)/10);
    deltam[atm3][atm2]=1-deltam[atm2][atm3];

	xi[0]=113;
    xi[1]=700;
    xi[2]=4;

    //load earnings profile from file profile.dta
	hcap = loadmat(indir+"profile.dta");
	if (isint(hcap)) oxwarning("file load failed "+indir+"profiles.dta");
    hcap=hcap[][1];
	
    }

 //given values, fill in vectors for At, L
exogenousfill(sc,pricefile,assetfile) {
  decl i,ii;

    //Initial Values
    CumC[]=213;
    Mb=590;


	cstate = reshape( <690,770,19180,0,.295>, D::P, 6);

	rho[]=1;

    for(i=1;i<D::P;++i) {
        phi[i]   = phi[i-1]   *(1+gammaphi*(1-deltaphi)^i);
	    Omegat[i]= Omegat[i-1]*(1+gammao*(1-deltao)^i);
	    theta[i] = theta[i-1] *(1+gammath*(1-deltath*deltath2^i)^i);
	   }
    Omega=Omegat;//set initial Omega to trend value - bring in damages in equil()
    decl j;
    pop[0][0]=88.5E6;
    for(j=1;j<D::G;++j) pop[0][j]=pop[0][j-1]./(1+gamman);
    for(i=1;i<D::P;++i) {
	   pop[i][0]=pop[i-1][0].*(1+gamman*(1-deltan)^(i));
       for(j=1;j<D::G;++j) pop[i][j]=pop[i-1][j-1];
       }
    N   = pop.*hcap';//N is units of effective labour supply per cohort per time period
    L   = sumr(N);  //L is effective labour supply in a time period
    N  *= 1E-6; //labour supply by cohort
    L  *= 1E-6; //aggregate labour supply.
    pop*= 1E-6;

    popshare=ones(D::P,D::G)./sumr(pop.*(1E6));  //each cohort's share of the total population
    //each agent in a particular cohort makes up popshare of the population
    //cohort share of pop is popshare.*pop.*10^6
    rentshare = permitshare = popshare;

	prices[][equity]=1*cumprod(constant(1,D::P,1));//start with a simple interest-free savings
    //prices[15:D::P-1][equity]=prices[15][equity];
    prices[][lab]=.004;
    prices[][res]=113.5;
    prices[][permits]=0;

    assets = loadmat(assetfile);
	prices = loadmat(pricefile);
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


 //procedure to assign parameters of the model for empirical calibration
parset(params) {
	gammao=params[0];
	deltao=params[1];
	Omegat[][]=params[2];

	gammath=params[3];
	deltath=params[4];
	theta[][]=params[5];

	gammaphi=params[6];
	deltaphi=params[7];
	phi[0][0]=params[8];

	exogenousfill(2,indir+"pricesfile.dat",indir+"assetsfile.dat"); //fill vectors of exogenous state variables and policy parameters

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
	println(" gammao=",params[0],"; //calibrated");
	println(" deltao=",params[1],"; //calibrated");
	println("Omegat[][]=",params[2],"; //calibrated");

	println(" ");
	println("gammaphi=",params[6],"; //calibrated");
	println("deltaphi=",params[7],"; //calibrated");
	println(" phi[0][0]=",params[8],"; //calibrated");
	println(" ");
	println(" ");
	println("gammath=",params[3],"; //calibrated");
	println("deltath=",params[4],"; //calibrated");
	println("theta[][]=",params[5],"; //calibrated");
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
	println(" gammao=",params[0],"; //calibrated");
	println(" deltao=",params[1],"; //calibrated");
	println("Omegat[][]=",params[2],"; //calibrated");

	println(" ");
	println("gammaphi=",params[6],"; //calibrated");
	println("deltaphi=",params[7],"; //calibrated");
	println(" phi[0][0]=",params[8],"; //calibrated");
	println(" ");
	println(" ");
	println("gammath=",params[3],"; //calibrated");
	println("deltath=",params[4],"; //calibrated");
	println("theta[][]=",params[5],"; //calibrated");
	println(" ");
	println("deltath2=1;");
	fclose("l");
	//limits used to keep parameters in possible ranges - positive TFP, theta, and emissions ratios, positive or negative initial growth rates, etc.
	}


calibprint() {
	decl i;
	for(i=1;i<D::P;++i)
		println(i,"  ",Omegat[i],"   ",Omega[i],"   ",cstate[i][],"     ",E[i],"     ",CumC[i],"     ",mextcost[i],"     ",E[i]*phi[i],"       ",L[i],"      ",double(sumr(pop[i][])),"      ",pop[i][0],"   ",K[i],"   ",U[i],"    ",GDP[i],"    ",cons[i][0],"   ",sumr(cons[i][].*pop[i][0]*(1E-6))[0][0],"   ",phi[i][0],"   ",theta[i][0],"   ",prices[i][0],"   ",prices[i][1],"   ",prices[i][2],"   ",prices[i][3],"   ",P[i][0],"   ",Pfirm[i][0],"   ",ctax[i][0],"   ",rentinc[i][0],"   ",permitinc[i],"   ",dividends[i],"   ",income[i],"   ",permit_transfers[i]);
		}


policy(permits_issued,carbon_tax,percap,sp,firm_allocation) {
	println("in policy ",permits_issued~carbon_tax~percap~sp~firm_allocation);

	P[:sp-1]= CarbPrice;//set initial permits to be non-binding
	if(permits_issued==CarbPrice) 
		//println("load dollar shadow price file");
		P[sp:] = loadmat(indir+"shadow_"+sprint(CarbPrice)+".dat")[sp:];
	else
		P[sp: ] = permits_issued;

	if(percap==1)
		P[sp:]=P[sp:]./sumr(pop[sp]).*sumr(pop[sp:]);  //per capita permits

	ctax[]=0;//the carbon tax rate is zero for early periods
	ctax[sp:]=carbon_tax;//set carbon tax for policy period
	Pfirm=P.*firm_allocation;
	}

		