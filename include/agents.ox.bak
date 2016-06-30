#include "humdim.oxh"

	
Utility(consumption) {	return sig1inv*pow(consumption,1-sigma) + Ubar;	}
	
MU(consumption) {	return 	pow(setbounds(consumption,.0001,+.Inf),-sigma);	}
	
Lifetime(consumption) {
	return beta.^range(0,rows(consumption)-1)*Utility(consumption); //CF Hardcoded .95 should have been beta??? 
	}
	
 //evaluate system of agent's euler equations
GV::agent_euler(sysval,startassets) {

	decl at = (EE|startassets|0.0);//first period endowment and death period zero added
	
	decl dt=(D|0).*at;//dividends per share times shareholdings

	at ~= (A[][1:]|0);
	
	decl
	  mu  = sumr(P.*at)[ : LS1];  //income
      mu -=  P[ : LS1][equity].*lag0(at[ : LS1][equity],-1)+dt[ : LS1];   //savings
      mu = MU(mu);

	sysval[0]   = -beta./P[ : LS1-1][equity];
	sysval[0] .*= (P[ 1 : LS1][equity]+D[ 1:LS1 ]);
	sysval[0] .*= mu[1:];
	sysval[0]  += mu[ : LS1-1]; 
	
	return !isnan(sysval[0]);  //CF this was not here before.	
	}
						
//equity distribution iteration
agents_problem(NTequil,itermax)	{
	decl fstart,firstper,maxper,masterp,topind,premax,padmax,curconv;
	decl born,startlife,Nconverge,euler; 	//CF removed decl of lifespan,
	decl storevec,assetvec,assetfill,consumption,consfill,utilityfill;
	decl np_age,sp_ext,step,MaxG;
	Nconverge=0;
	println("Beginning agent's problem");
	sp_ext=new matrix[D::G];//extra NG periods of share prices, holding last period ROR and dividend constant
	sp_ext[0]=shareprice[D::P-1]*prices[D::P-1][equity]-Inv::dividends[D::P-1];
	for(step=1;step<D::G;++step)  {
	  sp_ext[step]=sp_ext[step-1]*prices[D::P-1][equity]-Inv::dividends[D::P-1];
 	  }//end for loop													

	masterp = shareprice~(1E6*prices[][lab])~rentinc~permitinc;  //CF Avoids doing this repeatedly inside the loop
	MaxControl(itermax,0);  				//CF itermax was not being passed.

	MaxG = min(D::P,NTequil);
	
	for(born=-D::G+1; born<MaxG; ++born) {
		fstart=timer();
		if (!imod(born,10)) print(".");
		firstper = max(born,0);
		startlife = max(-born,0);
		GV::LS =  D::G-startlife;
		GV::LS1 = GV::LS-1;
		maxper  = firstper+GV::LS1;
		premax  = min(maxper,D::P-1);	
		topind = startlife+GV::LS1;
		np_age  = min(D::P-born,D::A)-1; //CF changed from		np_age= D::A;
		GV::P   = masterp[firstper:premax][];
		GV::D   = Inv::dividends[firstper:premax];//pass on relevant time period dividends per share
		GV::A   = (diagonal(assets[firstper:premax][startlife:np_age])')
						~hcap[startlife:np_age]
						~diagonal(rentshare[firstper:premax][startlife:np_age])'
						~diagonal(permitshare[firstper:premax][startlife:np_age])';		
		if (born==D::P-D::G-2) storevec=GV::A[][equity]/1000;   //CF original code stored all columns??
		padmax = born - (D::P-D::G) - 1;
		if(padmax>=0){
			GV::A |= 0~hcap[np_age+1:]~rentshare[D::P-1][np_age]~permitshare[D::P-1][np_age];
			GV::A[][equity]=storevec;		//try these as starting values - should be okay.
			GV::P   |=  sp_ext[:padmax]~reshape(masterp[D::P-1][1:],padmax+1,NPrices-1);
			GV::D   |=  constant(Inv::dividends[D::P-1],padmax+1,1);//pass on relevant time period dividends per share with fixed last period dividend through time
			}				
		if(GV::LS==1) {	//no choice
			cons[firstper][startlife]=sumr(GV::P[0][].*GV::A[0][]);
			}
		else  {  //CF: changed to else from "if(lifespan>1)"
			GV::EE = GV::A[0][equity];	//CF: removed sending lifespan in endow.
			assetvec=GV::A[1:GV::LS1][equity];  //CF replacment for line above

		GV::P  |= GV::P[rows(GV::P)-1][];
			curconv = SolveNLE(GV::agent_euler,&assetvec)>0;
			GV::agent_euler(&euler,assetvec);
			Nconverge += isnan(euler);
			GV::A[][equity]=((GV::EE|assetvec));
		GV::P  = GV::P[:rows(GV::P)-2][];
		GV::PxA = GV::P.*GV::A;
			consumption=sumr(GV::PxA)-GV::P[][equity].*lag0(GV::A[][equity],-1)
			            +(GV::D.*GV::A[][equity]);
			assetfill=assets[firstper:premax][startlife:];
			consfill=   cons[firstper:premax][startlife:];

			assetfill=setdiagonal(assetfill,endow[0]|submat(assetvec,0,np_age-1,0,0));
			consfill=setdiagonal(consfill,submat(consumption,0,np_age,0,0));

			assets[firstper:premax][startlife:]=assetfill;
			cons[firstper:premax][startlife:]=consfill;

		if(born>=0) {
			U[born]= Lifetime(consumption);
			permit_transfers[born][0]=GV::PxA[0][permits];
			income[born] = sumr(GV::PxA[0][]) + (GV::D[0].*GV::A[0][equity]);
			}
		}
	}//end of for loop over birth cohorts

	return Nconverge;
	}
