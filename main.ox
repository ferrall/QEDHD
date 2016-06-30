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

	policy(54,0.0485,0,48,0);
	scenario(2);
	equil("4");   //"2","2_50"

//	polisim();
//	calibration();
    }
