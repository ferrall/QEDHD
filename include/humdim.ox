﻿#include "setup.ox"
#include "firms.ox"
#include "agents.ox"

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

FirmEuler() {
	decl kold,xold,Kconv, Xconv,firm_crit, exstart, extstart, capstart, capchoice, extchoice;	
	extstart=shape(log(X),2*D::P,1);
	capstart=shape(K[1:],2*(D::P-1),1);//capital and shadow values

	Kconv = Xconv = FALSE;
	do {//solve extraction and production given rate of return
		kold=K;			xold=X;			E=X;
		newclimate(E);		//send sequence of emissions to temp routine
			//update climate damages
		Omega = Omegat./(d0+cstate[][stemp].*(b1 + b2*cstate[][stemp]));  //CF: simplified

		if (!Xconv) extchoice = extstart;
		Xconv = !SolveNLE(extract_nl,&extchoice,-1);			//println("Calculated Res Ext, given K with convergence crit ", SolveNLE(extract_nl,&extstart,-1));		
		X   = exp(extchoice[:D::P-1]);
		
		CumC=cumsum(lag0(X,1),1);

//		extrevenue=(theta.*Omega.*K.^alpha.*L.^(1-alpha-theta).*(phi.*X).^(theta));
//		extcost = (1E-3)*xi[0]*X+(1E-3)*(xi[1]*6000/(xi[2]+1))*(((CumC+X)/6000).^(xi[2]+1)-((CumC)/6000).^(xi[2]+1));

		if (!Kconv) capchoice = capstart;
		Kconv=!SolveNLE(invest_nl,&capchoice,-1);	//println("Calculated K, given res ext with convergence crit ",K_converge=SolveNLE(invest_nl,&capstart,-1));
		K[1:] = capchoice[:D::P-2] ;

		prices[][permits] = -extchoice[D::P : ];
		firm_crit = norm( (X|K)- (xold|kold),2);  				//CF almost same criterion firm_crit=meanc(meanr(((K~X)-(kold~xold))).^2); //mean square deviation
		if (firm_crit<firm_toler) println("Firm Converged"); 	//Convergence?, ",,"   ",X_converge,"   ",K_converge);
		} while(firm_crit>firm_toler);
	}

//main equilibrium procedure
equil(file_load,file_save) {
	decl crit_equil,agent_crit;
	decl adj,Pold;
	decl mpk,mpl,mpr;
	decl mcost;
	decl Aconv;
	decl disc = rplus1.^range(0,D::P-1)';  //Modified by CF 

	K = 19/rplus1*disc;					//initial guess of capital stock increases at the rate of discount
	X = 2*disc;							//initial guess of capital stock increases at the rate of discount

	// CF:  COMMENTED OUT THESE FILES TO AVOID POSSIBLE FAILURES
	//	K = loadmat(indir+"Kfile"~file_load~".dat");
	//	X=loadmat(indir+"Xfile"~file_load~".dat");

	assets=loadmat(indir+"assetsfile"~file_load~".dat");
	prices=loadmat(indir+"pricesfile"~file_load~".dat");

	Omega=Omegat;
	print_equil=0;
	shareprice= zeros(D::P,1);	//the price of one share in the last period in dollars
    do {	// loop over asset rate of return/shareholder discount factor.
		MaxControl(200,0);
		Pold=prices[][equity];
		FirmEuler();
		//		permitinc	     = (P!=Pfirm).*(X-Pfirm).*prices[][permits]*(1E12);  //CF:  Should !=  be .!= ????
		permitinc	     = (P .!= Pfirm).*(X-Pfirm).*prices[][permits]*(1E12);  //CF:  Should !=  be .!= ????
		firm_permit_value= -(Pfirm).*prices[][permits];
		//update share prices given dividends and rates of return
		shareprice[D::P-1]=1000;  //the price of one share in the last period in dollars
		for(timet=D::P-2;timet>=0;--timet)  
		  shareprice[timet]=(shareprice[timet+1]+dividends[timet+1])/prices[timet+1][equity];
		do {
			assets[1:][]=0;
     		Aconv = !agents_problem(1000);
            // println( Aconv ? "" : "Restart");
			} while(!Aconv);
		adj =spline( 1E-6 * sumr(assets.*pop) - 1 ,range(1,D::P)',0); //CF: simplified //locally smooth adjustment - keeps it from going pos/neg/pos
		agent_crit=norm(adj,2); 		//CF:  make criteria consisitent agent_crit=maxc(fabs(adj));
		prices[][equity] ./= (1+adj/50);  // Same as $50 // adjust rates of return
//		println(adj');
/*
		savemat(outdir+"xfile"~file_save~".dat",X);
		savemat(outdir+"kfile"~file_save~".dat",K);
		savemat(outdir+"assetsfile"~file_save~".dat",assets);
		savemat(outdir+"pricesfile"~file_save~".dat",prices);
*/
		crit_equil=norm((prices[][equity]-Pold)[:D::P-3],2); //firm crit will be small due to while loop above
		println("equil: ",crit_equil," ",agent_crit); //,"%c",{"old","new"},Pold~prices[][equity]);				// CF: println("equil_crit and conv flags",crit_equil,"  X=",X_converge,"   K=",K_converge,"  agent=",`onverge);
		} while (crit_equil+agent_crit>equil_toler);

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
	////savemat("outdir+paramsfile.dat",params);

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
	savemat(outdir+"xfile.dat",X);
	savemat(outdir+"kfile.dat",K);
	savemat(outdir+"assetsfile.dat",assets);
	savemat(outdir+"pricesfile.dat",prices);

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
	equil("4");   //"2","2_50"
	savemat(outdir+"shadow_"+sprint(CarbPrice)+".dat",X);

	println(ctax~prices[][permits]~X);
		
	//CF ?? This seems to be a duplicate.  Will run on first foreach()?
	policy(5.4,0,0,48,0);
	scenario(2);
	equil();
	savemat(outdir+"shadow_"+sprint(CarbPrice)+".dat",X);

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
