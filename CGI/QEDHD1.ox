#include "oxstd.h"
#import "solvenle"
#import <maximize>
#import "CGI"


decl beta, sigma, Iwealth, T; //number of periods
const decl dividend = 0.1;
decl cons,price,parlabs,parvals;

setparam(nm,par,defv) {
    decl loc= strfind(parlabs,nm);
    par[0] = loc>-1 ? parvals[loc] : defv;
    }

setparameters() {
    CGI::Initialize();
	println("\n\n----------------------------------------");
    [parlabs,parvals] = CGI::Parse();
    println("%c",{"Value"},"%r",parlabs,parvals);
	println("----------------------------------------");	
    CGI::Finalize();
    setparam("beta",&beta,0.96);
    setparam("sigma",&sigma,1.223);
    setparam("Iwealth",&Iwealth,3.0);
    setparam("T",&T,20);
    }

agent_euler(sysval, solution) {
	decl i,yt=matrix(Iwealth*1.04),yt1=<>;
	decl asset = solution[:T-2];
	 price=1+solution[T-1];   //hold asset returns constant over all t. (1XT)

	/*income vectors, yt=current period income and yt1=next period income.
    yt and yt1 are supposed to be (T-1)X1. I ignored labour income and resource firm dividends. */
	
	for(i=0;i<T-2;i++) {
	yt |= asset[i]*(price+dividend);
	}
	
	for(i=0;i<T-1;i++) {
	yt1 |= asset[i]*(price+dividend);
	}															

	//Euler equation. We have a system of T-1 equations. Remember that we had 1 equation when T=2.
	cons = yt1-(price*(asset[1:]|0));
	sysval[0] = -beta*((price+dividend)/price) *(cons).^-sigma; // (asset[1:]|0) since no bequest motive by the agent.
	sysval[0] = sysval[0]+(yt-(price*asset[:T-2])).^-sigma;
	sysval[0] |= sumc(asset) - 100;

	//println("asset return ", a_return,"dividends ", dividends, "asset ", asset, "income at t ", yt, "income at t+1 ", yt1);
	
	return TRUE;
	}
	
//Intro paragraph	
header(){
}

main () {
	// Asset vector is T-1 dimensions because we know the last stage has 0 demand
	// since the agent has no bequest motive. Also, T time stages yields T-1 steps.

    setparameters();		
	decl assetvec = constant(0.001,T-1,1), conval, euler;
	decl return_rate = 1;
	decl solution_vec = assetvec | return_rate;

	MaxControl(100000,-1);
	conval = SolveNLE(agent_euler, &solution_vec);	  //solving for vector of optimal asset holding.
	agent_euler(&euler, solution_vec);
	println("\nEuler: ", euler, "\nConvergence value = ", conval, "\n");
	println("\nOptimal asset demand and consumption: ", solution_vec[:T-2]~cons);
	println("\nSum of demands: ", sumc(solution_vec[:T-2]));
	println("\nReturn Rate: ", ((dividend/price)), "\nPrice ", price);
}
