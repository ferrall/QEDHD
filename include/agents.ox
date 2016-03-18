#include "humdim.oxh"

Utility(consumption) {
	return sig1inv*pow(consumption,1-sigma) + Ubar;
	}
	
MU(consumption) {
	return 	pow(setbounds(consumption,.0001,+.Inf),-sigma);
	}
	
Lifetime(consumption) {
	return beta.^range(0,rows(consumption)-1)*Utility(consumption); //CF Hardcoded .95 should have been beta??? 
	}
	
 //evaluate system of agent's euler equations
agent_euler(sysval,startassets) {
	//enum{equity,lab,res,permits}
	decl assetstemp=((endow[equity]|startassets)~assetspass[][lab:permits])|0.0;//first period endowment and death period zero added
	decl pricestemp=pricespass|pricespass[rows(pricespass)-1][];
	decl divstemp=(dividendspass|0).*assetstemp[][equity];//dividends per share times shareholdings
	
	decl
	  mu  = sumr(pricestemp.*assetstemp)[:lifespan-1];  //income
    mu -=  pricestemp[:lifespan-1][equity].*lag0(assetstemp[:lifespan-1][equity],-1)+divstemp[:lifespan-1];   //savings
    mu = MU(mu);

	sysval[0]   = -beta./pricestemp[:lifespan-2][equity];
	sysval[0] .*= (pricestemp[1:lifespan-1][equity]+dividendspass[1:lifespan-1]);
	sysval[0] .*= mu[1:];
	sysval[0]  += mu[:lifespan-2]; 
	
	sysval[0]*=10;//scale it for solution tightness
	//println("relevant measures - consumption, income, investment, dividends, euler",consumption~(sumr(pricestemp.*assetstemp))~(pricestemp[][equity].*lag0(assetstemp[][0],-1))~divstemp[][0]~(0|sysval[0]|0));
	//println(startassets~sysval[0]);
	return !isnan(sysval[0]);  //CF this was not here before.	
	}
						

//equity distribution iteration
agents_problem(itermax)	{
	decl fstart,firstper,maxper,masterp,topind,premax,padmax;
	decl born,startlife,Nconverge,euler; 	//CF removed decl of lifespan,
	decl storevec,assetvec,assetfill,consumption,consfill,utilityfill;
	decl np_age,sp_ext,step;
	Nconverge=0;
	println("Beginning agent's problem");
	sp_ext=new matrix[D::G];//extra NG periods of share prices, holding last period ROR and dividend constant
	sp_ext[0]=shareprice[D::P-1]*prices[D::P-1][equity]-dividends[D::P-1];
	for(step=1;step<D::G;++step)  {
	  sp_ext[step]=sp_ext[step-1]*prices[D::P-1][equity]-dividends[D::P-1];
 	  }//end for loop													

	masterp = shareprice~(1E6*prices[][lab])~rentinc~permitinc;  //CF Avoids doing this repeatedly inside the loop
	MaxControl(itermax,0);  //CF itermax was not being passed.
	for(born=-D::G+1;born<D::P;++born) {
		fstart=timer();
		if (!imod(born,10)) print(".");
		//CF Removed if() statement
		firstper = max(born,0);
		startlife = max(-born,0);
		lifespan =  D::G-startlife;
		maxper  = firstper+lifespan-1;
		premax  = min(maxper,D::P-1);	
		topind = startlife+lifespan-1;
		np_age  = min(D::P-born,D::A)-1; //CF changed from		np_age= D::A;
		pricespass   = masterp[firstper:premax][];
		dividendspass=dividends[firstper:premax];//pass on relevant time period dividends per share
		assetspass=  (diagonal(assets[firstper:premax][startlife:np_age])')
						~hcap[startlife:np_age]
						~diagonal(rentshare[firstper:premax][startlife:np_age])'
						~diagonal(permitshare[firstper:premax][startlife:np_age])';		
		if(born==D::P-D::G-2) storevec=assetspass[][equity]/1000;   //CF original code stored all columns??
		padmax = born - (D::P-D::G)-1;
		if(padmax>=0){
			assetspass |= 0~hcap[np_age+1:]~rentshare[D::P-1][np_age]~permitshare[D::P-1][np_age];
			assetspass[][equity]=storevec;		//try these as starting values - should be okay.
			pricespass   |=  sp_ext[:padmax]~reshape(masterp[D::P-1][1:],padmax+1,NPrices-1);
			dividendspass|=  constant(dividends[D::P-1],padmax+1,1);//pass on relevant time period dividends per share with fixed last period dividend through time
			}				
		consumption=sumr(pricespass.*assetspass);
		// assetvec=vec(assetspass[0][equity]);
		if(lifespan==1) {	//no choice
			cons[firstper][startlife]=consumption[0];
//			utils[firstper][startlife]=Utility(cons[firstper][startlife]);
			}
		else  {  //CF: changed to else from "if(lifespan>1)"
			endow=assetspass[0][];	//CF: removed sending lifespan in endow.
			assetvec=assetspass[1:lifespan-1][equity];  //CF replacment for line above
			Nconverge += SolveNLE(agent_euler,&assetvec);
			agent_euler(&euler,assetvec);
			assetspass[][equity]=((endow[equity]|assetvec));
			consumption=sumr(pricespass.*assetspass)-pricespass[][equity].*lag0(assetspass[][equity],-1)+(dividendspass.*assetspass[][equity]);
			assetfill=assets[firstper:premax][startlife:];
			consfill=   cons[firstper:premax][startlife:];
			assetfill=setdiagonal(assetfill,endow[0]|submat(assetvec,0,np_age-1,0,0));
			consfill=setdiagonal(consfill,submat(consumption,0,np_age,0,0));

			assets[firstper:premax][startlife:]=assetfill;
			cons[firstper:premax][startlife:]=consfill;

		if(born>=0) {
			U[born]= Lifetime(consumption);
			permit_transfers[born][0]=sumr(pricespass[0][3].*assetspass[0][3]);
			income[born][0]=sumr(pricespass[0][].*assetspass[0][])+(dividendspass[0][].*assetspass[0][0]);
			}
		}
	}//end of for loop over birth cohorts

//	utils = Utility(cons);

	println("Sum of Convergence flags from agent's problem (must be zero):",Nconverge);
	return Nconverge;
	}
