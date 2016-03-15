#include "humdim.oxh"

//procedure to set up variables
initialize() {
  xi=new matrix[3];

  prices=new matrix[D::P][4];

  E=  phi=  theta=
 income= permit_transfers= ctax=CumC= rho= GDP = Omegat=U= K= X= P= zeros(D::P,1);

  
  deltam=new matrix[Ntrans][Ntrans];

  N= pop=  assets=  cons=  utils=  newassets=
  popshare=  //share of population represented by each agent in each cohort
  rentshare= //share of resource rents
  permitshare= //share of permit rents
  	zeros(D::P,D::G);

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
    decl db,hcbase;
    db = new Database();
    db.LoadDta(ddir+"profile.dta",1,1,1);
    hcbase=db.GetAll();
    hcap=hcbase[][1];
	
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

    decl dbase;
    dbase = new Database();
    dbase.Load(assetfile);
    assets=dbase.GetAll();
    dbase = new Database();
    dbase.Load(pricefile);
    prices=dbase.GetAll();
    }

policy(permits_issued,carbon_tax,percap,sp,firm_allocation) {
	println("in policy ",permits_issued~carbon_tax~percap~sp~firm_allocation);

	P[:sp-1]= CarbPrice;//set initial permits to be non-binding
	P[sp: ] = permits_issued;
	if(permits_issued==CarbPrice) {
		//println("load dollar shadow price file");
		decl dbase = new Database();
		dbase.Load(ddir+"shadow_"+sprint(CarbPrice)+".dat");
		P[sp:]=dbase.GetAll()[sp:];
		}	

	if(percap==1)
		P[sp:]=P[sp:]./sumr(pop[sp]).*sumr(pop[sp:]);  //per capita permits

	ctax[]=0;//the carbon tax rate is zero for early periods
	ctax[sp:]=carbon_tax;//set carbon tax for policy period
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
	for(i=1;i<D::P;++i)
		println(i,"  ",Omegat[i],"   ",Omega[i],"   ",cstate[i][],"     ",E[i],"     ",CumC[i],"     ",mextcost[i],"     ",E[i]*phi[i],"       ",L[i],"      ",double(sumr(pop[i][])),"      ",pop[i][0],"   ",K[i],"   ",U[i],"    ",GDP[i],"    ",cons[i][0],"   ",sumr(cons[i][].*pop[i][0]*(1E-6))[0][0],"   ",phi[i][0],"   ",theta[i][0],"   ",prices[i][0],"   ",prices[i][1],"   ",prices[i][2],"   ",prices[i][3],"   ",P[i][0],"   ",Pfirm[i][0],"   ",ctax[i][0],"   ",rentinc[i][0],"   ",permitinc[i],"   ",dividends[i],"   ",income[i],"   ",permit_transfers[i]);
		}

 //evaluate system of agent's euler equations
invest_nl(sysval,capital) {
	decl euler,surplus,kt;
	decl sent=shape(capital,D::P-1,2);
	decl capchoice=(K[0]~0)|sent;   //CF simplified

	decl mpk=alpha.*Omega.*capchoice[][0].^(alpha-1).*L.^(1-alpha-theta).*(phi.*X).^(theta)+ones(D::P,1)-capchoice[][1];

	decl discount=beta.*ones(D::P,1)./prices[][equity];

	GDP=Omega.*capchoice[][0].^(alpha).*L.^(1-alpha-theta).*(phi.*X).^(theta);
	prices[][lab]=(1-alpha-theta).*GDP./L;
	prices[][res]=(theta).*GDP./X;

	surplus=Omega.*capchoice[][0].^(alpha).*L.^(1-alpha-theta).*(phi.*X).^(theta)-prices[][lab].*L-prices[][res].*X+prices[][permits].*Pfirm+capchoice[][0];
	dividends=surplus-lag0(capchoice[][0],-1)/(1-delta);//dividends per share from production in dollars
	dividends[D::P-1] -= capchoice[][0][D::P-1][0]./(1-delta); //force firm to maintain capital stock.

	//use .1 here - as long as it is not binding in equil it is fine...using exactly zero admits small neg values.
	kt=(dividends[:D::P-2].<=0.1).*(dividends[:D::P-2]-0.1)+(dividends[:D::P-2].>0.1).*capchoice[1:][1]+(capchoice[1:][1].>0).*capchoice[1:][1];
	euler=discount[:D::P-2].*mpk[1:][]-ones(D::P-1,1)./(1-delta);

	if(print_equil==1)
		println((0|euler)~(dividends)~capchoice~mpk~X);
	sysval[0]=(euler|kt);
	}

 //evaluate system of agent's euler equations
extract_nl(sysval,sent_vals) {
	decl i,mpr,revenue;
	decl sent=shape(sent_vals,D::P,2);
	decl extraction=exp(sent[][0]);			//constrain to be positive
	decl Cumtemp=cumsum(lag0(extraction,1),1);
	decl Cum_eps=cumsum(lag0(extraction+(1E-3|zeros(D::P-1,1)),1),1);
					  
	decl costs=(1E-3)*xi[0]*extraction+(1E-3)*(xi[1]*6000/(xi[2]+1))*(((Cumtemp+extraction)/6000).^(xi[2]+1)-((Cumtemp)/6000).^(xi[2]+1));
	//calculate the change in extraction costs for a small increase in extraction
	decl costs_eps=(1E-3)*xi[0]*extraction+(1E-3)*(xi[1]*6000/(xi[2]+1))*(((Cum_eps+extraction)/6000).^(xi[2]+1)-((Cum_eps)/6000).^(xi[2]+1));

	decl costinc=reversec(cumsum(reversec(costs_eps-costs),ones(D::P,1)*beta./reversec(prices[][equity])))./1E-3;

	decl opp_cost=costinc;
	//opp_cost is the current value of all future cost increases caused by a one unit change in extraction in that period

	decl marg_rev=(theta.*Omega.*K.^alpha.*L.^(1-alpha-theta).*(phi.*extraction).^(theta)./extraction)+1/Z*extraction.*(theta.*(theta-1).*Omega.*K.^alpha.*L.^(1-alpha-theta).*(phi.*extraction).^(theta)./(extraction.^2));//mkt power - prices are endogenous
	marg_rev+=sent[][1]-ctax;
	mextcost=(1E-3)*(xi[0]+xi[1].*((Cumtemp[][0]+extraction)/6000).^xi[2]);//marginal extraction cost

	decl euler=marg_rev-mextcost-1/Z*opp_cost;

	decl kt=(extraction.>P).*(extraction-P).^2+(extraction.<P).*sent[][1].^2+(sent[][1].>0).*sent[][1].^2;
	prices[][permits]=-sent[][1];

	mpr=theta.*Omega.*K.^(alpha).*L.^(1-alpha-theta).*(phi.*X).^(theta)./X;
	prices[][res]=mpr-prices[][permits]-ctax;
	resmarkup=mpr-prices[][permits]-ctax-mextcost;
	revenue=prices[][res].*extraction;
	rentinc=(revenue-costs)*(1E12);//total rent income in dollars to pass to agents' problem
	sysval[0]=euler|kt;
	if(print_equil==1)
		println("resource sector",extraction~prices[][res]~P~prices[][permits]~rentinc);
	}

 //evaluate system of agent's euler equations
new_assetsnl(sysval,startassets) {
	//enum{equity,lab,res,permits}
	decl assetstemp=((endow[equity]|startassets)~assetspass[][lab:permits])|0.0;//first period endowment and death period zero added
	decl pricestemp=pricespass|pricespass[rows(pricespass)-1][];
	decl divstemp=(dividendspass|0).*assetstemp[][equity];//dividends per share times shareholdings
	
	decl
	  utility  = sumr(pricestemp.*assetstemp)[:lifespan-1];  //income
	  utility -=  pricestemp[:lifespan-1][equity].*lag0(assetstemp[:lifespan-1][equity],-1)+divstemp[:lifespan-1];   //savings
	  utility = pow(setbounds(utility,.0001,+.Inf),-sigma);

	sysval[0]   = -beta./pricestemp[:lifespan-2][equity];
	sysval[0] .*= (pricestemp[1:lifespan-1][equity]+dividendspass[1:lifespan-1]);
	sysval[0] .*= utility[1:];
	sysval[0]  += utility[:lifespan-2]; 
	
	sysval[0]*=10;//scale it for solution tightness
	//println("relevant measures - consumption, income, investment, dividends, euler",consumption~(sumr(pricestemp.*assetstemp))~(pricestemp[][equity].*lag0(assetstemp[][0],-1))~divstemp[][0]~(0|sysval[0]|0));
	//println(startassets~sysval[0]);	
	return !isnan(sysval[0]);  //CF this was not here before.	
	}
						
//equity distribution iteration
agents_problem(itermax)	{
	decl fstart,firstper,maxper,lastper,masterp,topind;
	decl born,startlife,Nconverge,euler; 	//CF removed decl of lifespan,
	decl storevec,assetvec,assetfill,consumption,consfill,utilityfill;
	decl np_age,sp_ext,step;
	Nconverge=0;
	//println("Beginning agent's problem");
	sp_ext=new matrix[D::G];//extra NG periods of share prices, holding last period ROR and dividend constant
	sp_ext[0]=shareprice[D::P-1]*prices[D::P-1][equity]-dividends[D::P-1];
	for(step=1;step<D::G;++step)  {
	  sp_ext[step]=sp_ext[step-1]*prices[D::P-1][equity]-dividends[D::P-1];
 	  }//end for loop													

	//println(shareprice|sp_ext);

	masterp = shareprice~(1E6*prices[][lab])~rentinc~permitinc;
	
	for(born=-D::G+1;born<D::P;++born) {
		fstart=timer();
		//CF Removed if() statement
		firstper = max(born,0);
		startlife = max(-born,0);
		lifespan =  D::G-startlife;
		maxper  = firstper+lifespan-1;
		lastper = min(maxper,D::P-1);
		topind = startlife+lifespan-1;
		np_age= D::A;
		pricespass   = masterp[firstper:lastper][];
		dividendspass=dividends[firstper:lastper];//pass on relevant time period dividends per share
		assetspass=  (diagonal(assets[firstper:D::P-1][startlife:D::G-1])')
						~hcap[startlife:D::G-1]
						~diagonal(rentshare[firstper:D::P-1][startlife:D::G-1])'
						~diagonal(permitshare[firstper:D::P-1][startlife:D::G-1])';		
//		if(born<=D::P-D::G)	{			//prices paid on assets are equity price+per share dividend, wage rate, total res rents, total permit rents
			if(born==D::P-D::G-2) storevec=assetspass[][equity]/1000;   //CF original code stored all columns??
//			}
		if(born>D::P-D::G){
			np_age=D::P-born;
			assetspass |= 0~hcap[]~rentshare[D::P-1][np_age-1]~permitshare[D::P-1];
//			assetspass=diagonal(assets[firstper:D::P-1][startlife:D::G-1]' |zeros(D::G-np_age,1))
//			          ~hcap[startlife:topind]
//					  ~(diagonal(rentshare[firstper:D::P-1][startlife:np_age-1])'  |constant(rentshare[D::P-1][np_age-1],D::G-np_age,1)))
//					  ~(diagonal(permitshare[firstper:D::P-1][startlife:np_age-1])'|cosntant(permitshare[D::P-1],D::G-np_age,1) ) ;
			assetspass[][equity]=storevec;		//try these as starting values - should be okay.
			pricespass   |=  sp_ext[:D::G-np_age-1]~reshape(masterp[D::P-1][1:],D::G-np_age,NPrices-1);
			dividendspass|=  constant(dividends[D::P-1],D::G-np_age,1);//pass on relevant time period dividends per share with fixed last period dividend through time
			}				
		println("Born in period ",born," with lifespan ",lifespan," age at D::P=",np_age," padded with ",D::G-np_age," zeros. asset,price lengths",rows(pricespass)," ",rows(assetspass));
		//		println("assetpass ",assetspass,"price pass ",pricespass,"div ", dividendspass);
		//assets are equity, human capital, res rent allocation, and permit rent allocation
		consumption=sumr(pricespass.*assetspass);
		assetvec=vec(assetspass[0][0]);
		if(lifespan==1) {
			cons[firstper][startlife]=consumption[0];
			utils[firstper][startlife]=cons[firstper][startlife]^(1-sigma)/(1-sigma)+Ubar;
			//println("Consumption ",(cumsum(ones(rows(consumption),1),1)-1)~(sumr(pricespass.*assetspass))~consumption);
			}
		else  {  //CF: changed to else from "if(lifespan>1)"
			//  endow=lifespan~assetspass[0][];//row vector of all three endowments
			endow=assetspass[0][];	//CF: removed sending lifespan in endow.
			// CF:vec() redundant?			assetvec=vec(assetspass[1:lifespan-1][equity]);//vector of capital holdings, in partial shares
			assetvec=assetspass[1:lifespan-1][equity];  //CF replacment for line above
			MaxControl(1000,0);
			//CF:not needed?			euler=assetvec;
			Nconverge += SolveNLE(new_assetsnl,&assetvec);
			new_assetsnl(&euler,assetvec);
			assetspass[][equity]=((endow[equity]|assetvec));
			consumption=sumr(pricespass.*assetspass)-pricespass[][equity].*lag0(assetspass[][equity],-1)+(dividendspass.*assetspass[][equity]);
			assetfill=assets[firstper:lastper][startlife:];
			consfill=   cons[firstper:lastper][startlife:];
			assetfill=setdiagonal(assetfill,endow[0]|submat(assetvec,0,np_age-2,0,0));
			consfill=setdiagonal(consfill,submat(consumption,0,np_age-1,0,0));

			assets[firstper:lastper][startlife:]=assetfill;
			cons[firstper:lastper][startlife:]=consfill;

		if(born>=0) {
			U[born]=sumc((consumption.^(1-sigma)./(1-sigma)+Ubar).*.95.^(cumsum(ones(rows(consumption),1),1)-1));
			permit_transfers[born][0]=sumr(pricespass[0][3].*assetspass[0][3]);
			income[born][0]=sumr(pricespass[0][].*assetspass[0][])+(dividendspass[0][].*assetspass[0][0]);
			}
		}
	}//end of for loop over birth cohorts

	utils=cons.^(1-sigma)./(1-sigma)+Ubar;

	println("Sum of Convergence flags from agent's problem (must be zero):",Nconverge);
	return Nconverge;
	}

//évolution du modèle climatique étant donné les émissions.
newclimate(emissions)	{
	//use global D::Px6 matrix of climate values and receive a D::Px1 vector of emissions levels
	//fill a new D::Px6 matrix of updated climate values.
	decl i;
	//cstate-=<596.4,705,19200,0,0,0>;//convert to deviations from P.I. normals
	for(i=1;i<D::P;++i)	{
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
	K=19/1.005*cumprod(constant(1.005,D::P,1));//initial guess of capital stock increases at the rate of discount
	X=2*cumprod(constant(1.005,D::P,1));//initial guess of capital stock increases at the rate of discount
	decl ext_euler=X;
	decl inv_euler=K;

	dbase = new Database();
	dbase.Load(ddir+"Kfile"~file_load~".dat");
	K=dbase.GetAll();
	delete dbase;  //CF Andrew's code did not clean up.

	dbase = new Database();
	dbase.Load(ddir+"Xfile"~file_load~".dat");
	X=dbase.GetAll();
	delete dbase;  //CF Andrew's code did not clean up.

	dbase = new Database();
	dbase.Load(ddir+"assetsfile"~file_load~".dat");
	assets=dbase.GetAll();
	delete dbase;  //CF Andrew's code did not clean up.

	dbase = new Database();
	dbase.Load(ddir+"pricesfile"~file_load~".dat");
	prices=dbase.GetAll();
	delete dbase;  //CF Andrew's code did not clean up.

	Omega=Omegat;
    do {
		MaxControl(200,0);
		ror_old=prices[][equity];
		do {//solve extraction and production given rate of return
			kold=K;
			xold=X;
			E=X;
			newclimate(E);//send sequence of emissions to temp routine
			//update climate damages
			Omega=Omegat./(constant(d0,D::P,1).*(ones(D::P,1)+b1*cstate[][stemp]+b2*cstate[][stemp].^2));
			extstart=log(X)|zeros(D::P,1);
			//println("Calculated Res Ext, given K with convergence crit ", SolveNLE(extract_nl,&extstart,-1));
        	X_converge=SolveNLE(extract_nl,&extstart,-1);
			print_equil=0;
			extract_nl(&ext_euler,extstart);
			print_equil=0;
			extchoice=shape(extstart,D::P,2);
			X=exp(extchoice[][0]);
			CumC=cumsum(lag0(X,1),1);
			extrevenue=(theta.*Omega.*K.^alpha.*L.^(1-alpha-theta).*(phi.*X).^(theta));
			extcost=(1E-3)*xi[0]*X+(1E-3)*(xi[1]*6000/(xi[2]+1))*(((CumC+X)/6000).^(xi[2]+1)-((CumC)/6000).^(xi[2]+1));

			capstart=(K[1:D::P-1][])|zeros(D::P-1,1);//capital and shadow values
			//println("Calculated K, given res ext with convergence crit ",K_converge=SolveNLE(invest_nl,&capstart,-1));
			K_converge=SolveNLE(invest_nl,&capstart,-1);
			print_equil=0;
			invest_nl(&inv_euler,capstart);
			print_equil=0;
			capchoice=shape(capstart,D::P-1,2);
			K=(K[0][0]|(capchoice[][0]));
		
			//update omega
			//println("Updating climate given res ext");
			prices[][permits]=-extchoice[][1];
			firm_crit=meanc(meanr(((K~X)-(kold~xold))).^2); //mean square deviation
			xold=X;
			kold=K;
			println("Firm Convergence crit, X, K ",firm_crit[0][0],"   ",X_converge[0][0],"   ",K_converge[0][0]);
			} while(firm_crit>firm_toler);

		print_equil=0;
		extract_nl(&ext_euler,extstart);
		print_equil=0;
		permitinc=(P!=Pfirm).*(X-Pfirm).*prices[][permits]*(1E12);
		firm_permit_value=-(Pfirm).*prices[][permits];
		//update share prices given dividends and rates of return
		shareprice=ones(D::P,1)*1000;//the price of one share in the last period in dollars
		for(timet=D::P-2;timet>=0;--timet)  {
		  shareprice[timet][0]=(shareprice[timet+1][0]+dividends[timet+1])/prices[timet+1][equity];
	 	  }//end for loop
		do {
     		agent_converge=agents_problem(1000);
            if (agent_converge>0) {
                println("Unconverged assets - recalculating");
			     assets[1:][]=0;
                }
			} while(agent_converge>0);
		  adj=((sumr(assets.*pop)*(1E6))./(1E12))-ones(D::P,1);
		  adj[][]=spline(adj,cumsum(ones(D::P,1),1),0);//locally smooth adjustment - keeps it from going pos/neg/pos
		  agent_crit=maxc(fabs(adj[][]));
		  prices[][equity]./=(1+adj[][]/50);  // Same as $50 // adjust rates of return 
		  savemat(ddir+"xfile"~file_save~".dat",X);
		  savemat(ddir+"kfile"~file_save~".dat",K);
		  savemat(ddir+"assetsfile"~file_save~".dat",assets);
		  savemat(ddir+"pricesfile"~file_save~".dat",prices);
		  crit_equil=meanc(fabs((prices[][equity])-(ror_old))[:D::P-3][])+agent_crit; //firm crit will be small due to while loop above
		  println("equil_crit and conv flags",crit_equil,"  X=",X_converge,"   K=",K_converge,"  agent=",agent_converge);
		} while (crit_equil>equil_toler);

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

	exogenousfill(2,"pricesfile.dat","assetsfile.dat"); //fill vectors of exogenous state variables and policy parameters

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

calibration() {
	params=gammao|deltao|Omega[0][0]|gammath|deltath|theta[0][0];
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
	////savemat("ddir+paramsfile.dat",params);

	params=params|params_phi;

	println(" //Calibrated Parameters - Output");
	println(" ");
	println(" gammao=",params[0],"; //calibrated");
	println(" deltao=",params[1],"; //calibrated");
	println("Omegat[][]=",params[2],"; //calibrated");
	println(" ");
	println(" ");
	println("gammath=",params[3],"; //calibrated");
	println("deltath=",params[4],"; //calibrated");
	println("theta[][]=",params[5],"; //calibrated");
	println(" ");
	println("deltath2=1;");
	println(" ");
	println(" ");
	println("gammaphi=",params[6],"; //calibrated");
	println("deltaphi=",params[7],"; //calibrated");
	println(" phi[0][0]=",params[8],"; //calibrated");
	println(" ");
	println(" ");

}

//  Calibrated model output
caloutput() {


	policy(54,0,0,48,0);
	scenario(2);
	equil();
	fopen("calib_rev2.log","l");
	calibprint();
	fclose("l");
	savemat(ddir+"xfile.dat",X);
	savemat(ddir+"kfile.dat",K);
	savemat(ddir+"assetsfile.dat",assets);
	savemat(ddir+"pricesfile.dat",prices);

	}

//Policy Simulations
polisim() {	

	decl logfile,g;

	////policy(permits_issued,carbon_tax,percap,sp,firm_allocation)
	decl quotachoices=<540;54;500>,qc;
	decl percapchoices=<0;1>, pcc;
	decl firmchoices=<0;1>, fc;
	decl scchoices=<2>, sc;


	//find CarbPrice shadow value emissions
	policy(54,0.0485,0,48,0);
	scenario(2);
	equil();   //"2","2_50"
	savemat(ddir+"shadow_"+sprint(CarbPrice)+".dat",X);

	println(ctax~prices[][permits]~X);

	//CF ?? This seems to be a duplicate.  Will run on first foreach()?
	policy(5.4,0,0,48,0);
	scenario(2);
	equil();
	savemat(ddir+"shadow_"+sprint(CarbPrice)+".dat",X);

	println(ctax~prices[][permits]~X);

	foreach(qc in quotachoices[g])
		foreach( sc in scchoices )
			if(g) {  // CF ?? I think if - else not needed
				foreach(pcc in percapchoices)
					foreach( fc in firmchoices) {
						//println("permit allocation ",qc," per capita indication ",pcc," firm permits ",fc, " and climate scenario ",sc);
						policy(qc/10,0,pcc,48,fc);
						//println("permit allocation and firm permits",P~Pfirm);
						logfile="simuls_"~sprint(qc)~"_"~sprint(pcc)~"_"~sprint(fc)~"_"~sprint(sc)~".log";
						println("logfile is ",logfile);
						scenario(sc);
						println("t2x and temp change parameter ",t2x~t2);
						equil(sprint(sc),sprint(sc)~"_"~sprint(qc));
						fopen(logfile,"l");
						calibprint();
						fclose("l");
						}
				}
			else {
            //					//println("permit allocation ",qc," per capita indication ",pcc," firm permits ",fc, " and climate scenario ",sc);
				policy(54,0,0,48,0);
            //					//println("permit allocation and firm permits",P~Pfirm);
				logfile="simuls_"~sprint(qc)~"_"~sprint(0)~"_"~sprint(0)~"_"~sprint(sc)~".log";
				println("logfile is ",logfile);
				scenario(sc);
				println("t2x and temp change parameter ",t2x~t2);
				equil(sprint(sc),sprint(sc));
				fopen(logfile,"l");
				calibprint();
				fclose("l");
				}
			
		
		fopen("cohorts.log","l");
		println(cumsum(ones(D::P,1),1)~pop[][0]);
		fclose("l");
	}
