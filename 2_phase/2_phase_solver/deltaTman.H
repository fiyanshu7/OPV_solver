phiEqniRes = phiEqnSol.initialResidual();
//phiEqnfRes = phiEqnSol.finalResidual();

mu_eEqniRes = mu_eEqnSol.initialResidual();
//mu_eEqnfRes = mu_eEqnSol.finalResidual();

mu_hEqniRes = mu_hEqnSol.initialResidual();
//mu_hEqnfRes = mu_hEqnSol.finalResidual();

//Selecting the largest redisual to control deltaT
if (phiEqniRes>mu_eEqniRes)
    lim_res = phiEqniRes;
else
    lim_res = mu_eEqniRes;

if(lim_res < mu_hEqniRes)
    lim_res = mu_hEqniRes;

if(lim_res<5.0e-5 && t_gap%50==0)
{
     runTime.setDeltaT(1.1*runTime.deltaTValue());
     if (runTime.deltaTValue()>20.48)
         runTime.setDeltaT(20.48);
}


t_gap++;

// Info  << " phiEqninitResidual = " << phiEqniRes << endl;
// Info  << " phiEqnfinResidual  = " << phiEqnfRes << endl;
// Info  << " mu_eEqninitResidual = " << mu_eEqniRes << endl;
// Info  << " mu_eEqnfinResidual  = " << mu_eEqnfRes << endl;
// Info  << " mu_hEqninitResidual = " << mu_hEqniRes << endl;
// Info  << " mu_hEqnfinResidual  = " << mu_hEqnfRes << endl;
