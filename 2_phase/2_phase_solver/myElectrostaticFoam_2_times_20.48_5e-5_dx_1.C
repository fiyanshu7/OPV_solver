/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    electrostaticFoam

Description
    Solver for electrostatics.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    #include "fermiEqn.H"
	
    if(runTime.value()==0)
    {
	//runTime.setDeltaT(1e-1);
	Info<< nl << "deltaT = " << runTime.deltaTValue() << nl << endl;

    	scalar t=0;
	for(t=0;t<=10;t++)
	{         
        	solve
        	(
            	sig1*4.0*fvm::ddt(interface) - sig*4.0*fvm::laplacian(interface) + (0.25)*(2*interface + 4*interface*interface*interface - 6*interface*interface)
        	);
    	}
    }
    
    volScalarField comp_e_init(interface*c_e_P3HT +  (1.0 - interface)*c_e_PCBM);
    volScalarField comp_h_init(interface*c_h_P3HT +  (1.0 - interface)*c_h_PCBM);
    
    volScalarField ni_squared(comp_e_init*comp_h_init);
    background = -1.0*(q_e*comp_e_init + q_h*comp_h_init);
    
    if(runTime.value()==0)
    {
        Info<< nl << "Solving for mu" << endl;
        #include "muCalc.H"
    }

    scalar phiEqniRes;
    scalar mu_eEqniRes;
    scalar mu_hEqniRes;
    scalar lim_res;
 
    label t_gap=0;
    
    Info<< "\nStarting Time loop\n" << endl;
    runTime.setDeltaT(20.48);    
			 
    Info<< "Time = " << runTime.timeName() << nl << endl;	
  
    while (runTime.loop())
    {
        Info<< "deltaT = " << runTime.deltaTValue() << endl;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        comp_e = interface*exp((mu_e - E_e_P3HT)/(R*T)) + (1.0-interface)*exp((mu_e - E_e_PCBM)/(R*T));
        comp_h = interface*exp((mu_h - E_h_P3HT)/(R*T)) + (1.0-interface)*exp((mu_h - E_h_PCBM)/(R*T));
 
        rho = (q_e*comp_e + q_h*comp_h + background)*faraday;
         
        volScalarField M_e(interface*D_e_P3HT*exp((mu_e - E_e_P3HT)/(R*T))/(R*T)  +  (1 - interface)*D_e_PCBM*exp((mu_e - E_e_PCBM)/(R*T))/(R*T));     
        volScalarField M_h(interface*D_h_P3HT*exp((mu_h - E_h_P3HT)/(R*T))/(R*T)  +  (1 - interface)*D_h_PCBM*exp((mu_h - E_h_PCBM)/(R*T))/(R*T)); 
         
        volScalarField chi_e(interface*exp((mu_e - E_e_P3HT)/(R*T)) + (1.0 - interface)*exp((mu_e - E_e_PCBM)/(R*T)));  
        volScalarField chi_h(interface*exp((mu_h - E_h_P3HT)/(R*T)) + (1.0 - interface)*exp((mu_h - E_h_PCBM)/(R*T)));
         
        fvScalarMatrix phiEq
        (
            sig*fvm::laplacian(phi) == -1.0*(rho/epsilon0)
        );
        
        const solverPerformance phiEqnSol = phiEq.solve();
        
        /*
            //Another scheme for solving div.(M.grad(emu))
            solve
            (
                sig1*chi_e*fvm::ddt(emu_e) - sig*(fvm::laplacian(M_e, emu_e) - (fvc::grad(M_e) & fvc::grad(emu_e))) - 6*interface*(1-interface)*(G0-                                                R0*M_h*(comp_e*comp_h - ni_squared))
            );
        */
        
        fvScalarMatrix mu_eEq
        (
            sig1*chi_e*fvm::ddt(mu_e) == sig*(fvm::laplacian(M_e, mu_e) + (fvc::laplacian(M_e, q_e*faraday*phi))) + 6*interface*(1-interface)*G0 - R0*M_h*(comp_e*comp_h - ni_squared)
        );
        
        const solverPerformance mu_eEqnSol = mu_eEq.solve();
        
        
        fvScalarMatrix mu_hEq
        (
            sig1*chi_h*fvm::ddt(mu_h) == sig*(fvm::laplacian(M_h, mu_h) + (fvc::laplacian(M_h, q_h*faraday*phi))) + 6*interface*(1-interface)*G0 - R0*M_h*(comp_e*comp_h - ni_squared)
        );
        
        const solverPerformance mu_hEqnSol = mu_hEq.solve();
        
	/*			
        if(runTime.value()>=1.0 && runTime.deltaTValue()<20.48)
        {
            #include "deltaTman.H"
        }
	*/

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
