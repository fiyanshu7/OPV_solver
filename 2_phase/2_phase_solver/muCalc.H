dimensionedScalar error_1
(
"error_1",
dimensionSet(0,0,0,0,0,0,0),
scalar(0.0)
); 
dimensionedScalar error_2
(
"error_2",
dimensionSet(0,0,0,0,0,0,0),
scalar(0.0)
);
dimensionedScalar dcbydmu_1
(
"dcbydmu_1",
dimensionSet(0,0,0,0,0,0,0),
scalar(1.0)
);
dimensionedScalar dcbydmu_2
(
"dcbydmu_2",
dimensionSet(0,0,0,0,0,0,0),
scalar(1.0)
);
dimensionedScalar dmu_1
(
"dmu_1",
dimensionSet(0,0,0,0,0,0,0),
scalar(0.0)
);
dimensionedScalar dmu_2
(
"dmu_2",
dimensionSet(0,0,0,0,0,0,0),
scalar(0.0)
);

dimensionedScalar tolmu
(
"tolmu",
dimensionSet(0,0,0,0,0,0,0),
scalar(1e-10)
);
//  scalar error_1=0;
//  scalar error_2=0;
//  scalar dmu_1=0;
//  scalar dmu_2=0;
//  scalar dcbydmu_1=0;
//  scalar dcbydmu_2=0;
    scalar iter          =  0;
    scalar max_iter      =  500;


forAll(interface, cellI) 
{	
    //Ideal solution model formulation
	
    if (interface[cellI] == scalar(1.0))
    {
        mu_e[cellI] = E_e_P3HT.value() + R.value()*T.value()*Foam::log(comp_e_init[cellI]);
        mu_h[cellI] = E_h_P3HT.value() + R.value()*T.value()*Foam::log(comp_h_init[cellI]);
    }
    
    else if (interface[cellI] == scalar(0.0))
    {
       mu_e[cellI] = E_e_PCBM.value() + R.value()*T.value()*Foam::log(comp_e_init[cellI]);
       mu_h[cellI] = E_h_PCBM.value() + R.value()*T.value()*Foam::log(comp_h_init[cellI]);
    }
    
//     mu_e[cellI] = (interface[cellI])*(E_e_P3HT.value() + R.value()*T.value()*Foam::log(comp_e_init[cellI])) + (1.0-interface[cellI])*(E_e_PCBM.value() + R.value()*T.value()*Foam::log(comp_e_init[cellI]));
// 
//     mu_h[cellI] = (interface[cellI])*(E_h_P3HT.value() + R.value()*T.value()*Foam::log(comp_h_init[cellI])) + (1.0-interface[cellI])*(E_h_PCBM.value() + R.value()*T.value()*Foam::log(comp_h_init[cellI]));
	    
    else
    {
        mu_e[cellI] =  0.0;
        mu_h[cellI] =  0.0;

        iter = 0;
        
        while(iter < max_iter) 
        {
            //Newton-Raphson implementation ----f(x+h) = f(x) + hf'(x)------ Here, f(x+h)-f(x) = error;   f'(x)=dcdmu;
               
            error_1     = comp_e_init[cellI] - (Foam::exp((mu_e[cellI] - E_e_P3HT.value())/(R*T))*interface[cellI] + Foam::exp((mu_e[cellI] - E_e_PCBM.value())/(R*T))*(1.0-interface[cellI]));
            error_2     = comp_h_init[cellI] - (Foam::exp((mu_h[cellI] - E_h_P3HT.value())/(R*T))*interface[cellI] + Foam::exp((mu_h[cellI] - E_h_PCBM.value())/(R*T))*(1.0-interface[cellI]));
            
//          Info<< "error_1 = " << mag(error_1) << endl;
//          Info<< "error_2 = " << mag(error_2) << endl;
//          Info<< "error_1 = " << error_1.value() << endl;
//          Info<< "error_2 = " << error_2.value() << endl;
            
            if ((mag(error_1) < tolmu) && (mag(error_2) < tolmu))
	    {
                break;
	    }
            
            dcbydmu_1   = (1.0/(R*T))*(Foam::exp((mu_e[cellI] - E_e_P3HT.value())/(R*T))*interface[cellI] + Foam::exp((mu_e[cellI] - E_e_PCBM.value())/(R*T))*(1.0-interface[cellI])); 
            dmu_1       = error_1/dcbydmu_1;
            mu_e[cellI]     += dmu_1.value();
            
            dcbydmu_2   = (1.0/(R*T))*(Foam::exp((mu_h[cellI] - E_h_P3HT.value())/(R*T))*interface[cellI] + Foam::exp((mu_h[cellI] - E_h_PCBM.value())/(R*T))*(1.0-interface[cellI])); 
            dmu_2       = error_2/dcbydmu_2;
            mu_h[cellI] += dmu_2.value();
            
            iter++;
             
            if(iter==max_iter)
            {
                Pout<< "Solution not converging!!!"<< endl;
                exit(1);
            }
//             printf("iter=%ld\n",iter);
//             printf("error_1=%le\terror_2=%le\n",error_1,error_2);
//             getchar();
        }
    } 
}
