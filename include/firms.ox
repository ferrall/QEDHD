#include "humdim.oxh"

Output(K,L,X) {
	return Omega.*K.^(alpha) .* L.^(1-alpha-theta) .* (phi.*X).^theta; 
	}
Ecosts(E,T) {
	return C1m *( xi[0]*E+  C6K*xi[1]/(xi[2]+1)*( ((T+E)/C6K).^(xi[2]+1)-(T/C6K).^(xi[2]+1)) );
	}
Mcosts(ET) {
	return C1m*(xi[0]+xi[1]*(ET/C6K).^xi[2]);
	}

Inv::Production() {
	myGDP = Output(myK,L,X);
	mpk=  alpha*myGDP./myK + 1 - mySH;
	dividends  = alpha*myGDP + prices[][permits].*Pfirm + myK;  //CF 
	dividends  -= del1inv*lag0(myK,-1);		//dividends per share from production in dollars
	dividends[D::Pm1] -= del1inv*myK[D::Pm1]; //CF:deleted redundant [][0]???  force firm to maintain capital stock.
	}
	
 //evaluate system of agent's euler equations
Inv::nl(sysval,capital) {
	decl myT = idiv(rows(capital),2);
	myK = K[0] |capital[:myT-1] | (myT<D::Pm1 ? K[myT+1:] : <>);
	mySH =    0|capital[myT:]   | (myT<D::Pm1 ? SH[myT+1:]: <>); 

	Inv::Production();

	//CF version of code below.  Not sure I understand
	kt   = dividends[:myT-1]-kt_toler;
	kt   = (kt.<=0.0) .? kt .: mySH[1:myT];
	kt  += setbounds(mySH[1:myT],0.0,+.Inf);
	
	sysval[0] = beta*(mpk[1:myT]./prices[:myT-1][equity]) - del1inv;

	if (print_equil==1)  println((0|sysval[0])~(dividends)~myK~mySH~mpk~X);
	
	sysval[0] |= kt;
	
	return !isnan(sysval[0]);
	}

 //evaluate system of agent's euler equations
Ext::nl(sysval,sent_vals) {
	decl myT = idiv(rows(sent_vals),2);

	myX    = exp(sent_vals[:myT-1]) | (myT<D::P ? X[myT:] : <>);			//constrain to be positive
	xprice = sent_vals[myT:] | (myT<D::P ? -prices[myT:][permits] : <>);

	decl Cumtemp=cumsum(lag0(myX,1),1),
		 Cum_eps=cumsum(lag0(myX+(C1m|zeros(D::Pm1,1)),1),1),					  
		//calculate the change in extraction costs for a small increase in extraction
		costs= Ecosts(myX,Cumtemp),
		costs_eps=	Ecosts(myX,Cum_eps),
        opp_cost=reversec(cumsum(reversec(costs_eps-costs),beta./reversec(prices[][equity])))./ C1m;	//opp_cost is the current value of all future cost increases caused by a one unit change in extraction in that period

	myout = Output(K,L,myX);
	
	decl marg_rev = theta.*myout./myX.*(1 + (1/Z)*(theta-1) ) + xprice-ctax;

	mextcost = Mcosts(Cumtemp+myX); 
	
	sysval[0] = (marg_rev-mextcost-(1/Z)*opp_cost)[:myT-1];

	kt = (myX - P)[:myT-1];
	kt = (kt.>0.0) .? sqr(kt) .:  sqr(xprice[:myT-1]);	
	kt += sqr( setbounds(xprice[:myT-1],0.0,.Inf) );

	sysval[0] |= kt;
	
//	if(print_equil==1)	println("resource sector",extraction~prices[][res]~P~prices[][permits]~rentinc);
		
	return !isnan(sysval[0]);
	}

//évolution du modèle climatique étant donné les émissions.
newclimate(emissions)	{
	//use global D::Px6 matrix of climate values and receive a D::Px1 vector of emissions levels
	//fill a new D::Px6 matrix of updated climate values.
	decl i;
	//cstate-=<596.4,705,19200,0,0,0>;//convert to deviations from P.I. normals
	for(i=1;i<D::P;++i)	{
		cstate[i][atm1]=emissions[i-1]+deltam[atm1][atm1]*cstate[i-1][atm1]+deltam[atm2][atm1]*cstate[i-1][atm2];
		cstate[i][atm2]=               deltam[atm2][atm2]*cstate[i-1][atm2]+deltam[atm1][atm2]*cstate[i-1][atm1]+deltam[atm3][atm2]*cstate[i-1][atm3];
		cstate[i][atm3]=               deltam[atm3][atm3]*cstate[i-1][atm3]+deltam[atm2][atm3]*cstate[i-1][atm2];
		//cstate+=<596.4,705,19200,0,0,0>;  //transform back to levels
		//c_new=c_new;//
		cstate[i][stemp]=t4*cstate[i-1][otemp]+t1*cstate[i-1][stemp]+t2*log((cstate[i-1][atm1])/596.4);
		cstate[i][otemp]=cstate[i-1][otemp]+t3*(cstate[i-1][stemp]-cstate[i-1][otemp]);
		}
	}

FirmEuler(NT) {

	decl lagXK, Kconv, Xconv,firm_crit, extstart, capstart, capchoice, extchoice;	

	extstart = log(X[:NT-1]) | -prices[:NT-1][permits];
	
	if (NT>1) capstart = K[1:NT-1] | SH[1:NT-1]; 	//capital and shadow values

	Kconv = Xconv = FALSE;
	
	do {//solve extraction and production given rate of return
		lagXK = X|K;
		newclimate(E=X);		//send sequence of emissions to temp routine

		//update climate damages
		Omega = Omegat./(d0+cstate[][stemp].*(b1 + b2*cstate[][stemp]));  //CF: simplified

		if (!Xconv) extchoice = extstart;
		Xconv = !SolveNLE(Ext::nl,&extchoice,-1);			//println("Calculated Res Ext, given K with convergence crit ", SolveNLE(extract_nl,&extstart,-1));		

		X[:NT-1]   = exp(extchoice[:NT-1]);
		prices[:NT-1][permits] = -extchoice[NT: ];

		if (NT>1) {
			if (!Kconv) capchoice = capstart;
			Kconv=!SolveNLE(Inv::nl,&capchoice,-1);	//println("Calculated K, given res ext with convergence crit ",K_converge=SolveNLE(invest_nl,&capstart,-1));
			K[1:NT-1] = capchoice[:NT-2] ;
			SH[1:NT-1] = capchoice[NT-1:];
			}
			
		GDP = Output(K,L,X);
		prices[][res]= theta.*GDP./X -prices[][permits]-ctax;
		//		resmarkup=prices[][res]-mextcost;
		
		rentinc= prices[][res].*X - 1E12*Ecosts(X,cumsum(lag0(X,1),1)) ;//total rent income in dollars to pass to agents' problem		
		prices[][lab]= (1-alpha-theta).*GDP./L;

		firm_crit = norm( (X|K)- lagXK,2);  				//CF almost same criterion firm_crit=meanc(meanr(((K~X)-(kold~xold))).^2); //mean square deviation

		if (firm_crit<firm_toler) println("Firm Converged"); 	//Convergence?, ",,"   ",X_converge,"   ",K_converge);
		} while(firm_crit>firm_toler);
	}
	