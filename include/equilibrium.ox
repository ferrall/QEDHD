#include "humdim.oxh"

struct EQ {
	static decl Pold, Aold, Aconv,ExDemand;
	static Step(NT);
	}

EQ::Step(NT) {
	decl timet; 
	MaxControl(200,0);
	Pold=prices[:NT-1][equity];
	Aold = assets;
	FirmEuler(NT);		
	permitinc	     = (P .!= Pfirm).*(X-Pfirm).*prices[][permits]*(1E12);  //CF:  Should !=  be .!= ????
	firm_permit_value= -(Pfirm).*prices[][permits];
	for(timet=D::Pm2;timet>=0;--timet)  
		  shareprice[timet]=(shareprice[timet+1]+Inv::dividends[timet+1])/prices[timet+1][equity];	
	Aconv = !agents_problem(NT,1000);
	ExDemand = 1E-6*sumr(assets[:NT-1][].*pop[:NT-1][]) - 1;
	}

//main equilibrium procedure
equil(file_load,file_save) {
	decl crit_equil,agent_crit,nt,adj;
	decl mpk,mpl,mpr;
	decl mcost;
	decl Aconv;
	decl disc = rplus1.^range(0,D::Pm1)';  //Modified by CF 

	SH = K = 19/rplus1*disc;			//initial guess of capital stock increases at the rate of discount
	X = 2*disc;							//initial guess of capital stock increases at the rate of discount

	assets= loadmat(indir+"assetsfile"~file_load~".dat");
	prices= loadmat(indir+"pricesfile"~file_load~".dat");

	prices[][equity] = double(meanc(prices[][equity]));
	assets=reshape(meanc(assets),D::P,columns(assets));
	
	Omega=Omegat;
	print_equil=0;
	shareprice= zeros(D::Pm1,1)|1000;	//the price of one share in the last period in dollars

	do {	

		EQ::Step(D::P);

		for (nt=2;nt<=4;++nt) {
			do {
				EQ::Step(nt);
//				adj =spline( EQ::ExDemand ,range(1,nt)',0);
				agent_crit=norm(adj=EQ::ExDemand,2);
				prices[:nt-1][equity] ./= setbounds(1+adj,.85,1.05);
				println(EQ::Pold~prices[:nt-1][equity]~EQ::ExDemand);
				crit_equil=norm(prices[:nt-1][equity]-EQ::Pold,2);
				println(nt," ",agent_crit," ",crit_equil,adj');
				} while(FALSE); //agent_crit+ crit_equil> equil_toler
			}
			
         // println( Aconv ? "" : "Restart");
		//adj =spline( 1E-6 * sumr(assets.*pop) - 1 ,range(1,D::P)',0); //CF: simplified //locally smooth adjustment - keeps it from going pos/neg/pos
		//agent_crit=norm(adj,2); 		//CF:  make criteria consisitent agent_crit=maxc(fabs(adj));
		//prices[][equity] ./= (1+adj/50);  // Same as $50 // adjust rates of return

//		crit_equil=norm((prices[][equity]-Pold)[:D::P-3],2); //firm crit will be small due to while loop above

//		//println("equil: ",crit_equil," ",agent_crit); //,"%c",{"old","new"},Pold~prices[][equity]);				// CF: println("equil_crit and conv flags",crit_equil,"  X=",X_converge,"   K=",K_converge,"  agent=",`onverge);
//		println(Pold~prices[][equity]);

		} while (FALSE); //(crit_equil+agent_crit>equil_toler);
	}
