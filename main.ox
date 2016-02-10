#import "humdim"     // you have to include/ to the include path for this to work

main(){
    decl starttime,fintime,params,pol,sc;
    decl ptest,j,over,i;
    starttime=timer();
    NG=60;  //number of generations
    NP=220*(NG/60);//number of periods
    pstart=48;
    decl dbase = new Database();
    dbase.LoadObs("edata.txt",5,60,0,1,1,1);
    data=dbase.GetAll();
    initialize(); //set up all variables
    parameterize();//parameterize the model
    exogenousfill(2,"data/pricesfile.dat","data/assetsfile.dat"); //fill vectors of exogenous state variables and policy parameters
    Z=2;//set level of market power in resource sector to acheive approx \$100 per ton markup

    println(deltam);
    fintime=timer();
    println("Total Run Time ",timespan(starttime,fintime));
    }
