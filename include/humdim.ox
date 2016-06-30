#include "setup.ox"
#include "firms.ox"
#include "agents.ox"
#include "equilibrium.ox"

calibration() {
	params=gammao|deltao|Omega[0][0]|gammath|deltath|theta[0][0];
	decl params_phi=gammaphi|deltaphi|phi[0][0];
	limit_low=(params|params_phi) * .9;
	limit_low[3][0]=gammath*1.1;//gammath is negative
	limit_high=(params|params_phi) * 1.1;
	limit_high[3][0]=gammath * .9;//gammath is negative

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
