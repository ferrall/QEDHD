#include "humdim.oxh"

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
	return !isnan(sysval[0]);
	}

Ecosts(E,T) {
	return C1m *( xi[0]*E+  C6K*xi[1]/(xi[2]+1)*( ((T+E)/C6K).^(xi[2]+1)-(T/C6K).^(xi[2]+1)) );
	}

Mcosts(ET) {
	return C1m*(xi[0]+xi[1]*(ET/C6K).^xi[2]);
	}
	
 //evaluate system of agent's euler equations
extract_nl(sysval,sent_vals) {
	decl i,mpr,revenue;
	decl sent=shape(sent_vals,D::P,2);
	decl extraction=exp(sent[][0]);			//constrain to be positive
	decl Cumtemp=cumsum(lag0(extraction,1),1);
	decl Cum_eps=cumsum(lag0(extraction+(C1m|zeros(D::P-1,1)),1),1);
					  
	//calculate the change in extraction costs for a small increase in extraction
	decl costs= Ecosts(extraction,Cumtemp); 	//CF:  (1E-3)*xi[0]*extraction+(1E-3)*(xi[1]*6000/(xi[2]+1))*(((Cumtemp+extraction)/6000).^(xi[2]+1)-((Cumtemp)/6000).^(xi[2]+1));
	decl costs_eps=	Ecosts(extraction,Cum_eps);	// CF (1E-3)*xi[0]*extraction+(1E-3)*(xi[1]*6000/(xi[2]+1))*(((Cum_eps+extraction)/6000).^(xi[2]+1)-((Cum_eps)/6000).^(xi[2]+1));

	decl costinc=reversec(cumsum(reversec(costs_eps-costs),ones(D::P,1)*beta./reversec(prices[][equity])))./ C1m;

	decl opp_cost=costinc;
	//opp_cost is the current value of all future cost increases caused by a one unit change in extraction in that period

	decl marg_rev=(theta.*Omega.*K.^alpha.*L.^(1-alpha-theta) .* (phi.*extraction).^(theta)./extraction) + 1/Z*extraction.*(theta.*(theta-1).*Omega.*K.^alpha.*L.^(1-alpha-theta).*(phi.*extraction).^(theta)./(extraction.^2));//mkt power - prices are endogenous
	marg_rev+=sent[][1]-ctax;
	
	
	mextcost = Mcosts(Cumtemp+extraction); //	mextcost=C1m*(xi[0]+xi[1].*((Cumtemp[][0]+extraction)/6000).^xi[2]);//marginal extraction cost
	
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
	return !isnan(sysval[0]);
	}
