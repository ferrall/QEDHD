#import "include/humdim"    

main(){
    decl starttime,fintime,params,pol,sc;
    decl ptest,over;
    starttime=timer();
	//Not used??    pstart=48;
    decl dbase = new Database();
    dbase.LoadObs("edata.txt",5,60,0,1,1,1);
    data=dbase.GetAll();
    initialize(); //set up all variables
    parameterize();//parameterize the model
    exogenousfill(2,indir+"pricesfile.dat",indir+"assetsfile.dat"); //fill vectors of exogenous state variables and policy parameters
    Z=2;//set level of market power in resource sector to acheive approx \$100 per ton markup

	polisim();
//	calibration();
//    println(deltam);
//    fintime=timer();
//    println("Total Run Time ",timespan(starttime,fintime));
    }
