


volScalarField T0
(
    IOobject
    (
        "T0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField T1
(
    IOobject
    (
        "T1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField T2
(
    IOobject
    (
        "T2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Info<< "Reading field phi\n" << endl;
// volScalarField phi
// (
//     IOobject
//     (
//         "phi",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );
// 
// 
// volScalarField mu_e
// (
//     IOobject
//     (
//         "mu_e",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
//     //(E_e_P3HT+R*T*log(comp_e))*interface+(1-interface)*E_e_PCBM+R*T*log(comp_e)
// );
// 
// volScalarField mu_h
// (
//     IOobject
//     (
//         "mu_h",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
//     //(E_h_P3HT+R*T*log(comp_h))*interface+(1-interface)*E_h_PCBM+R*T*log(comp_h)
// );

//     volScalarField emu_e
//     (
//         IOobject
//         (
//             "emu_e",
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         //mesh
//        mu_e + q_e*faraday*phi
//     );
//     
//     volScalarField emu_h
//     (
//         IOobject
//         (
//             "emu_h",
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         //mesh
//         mu_h + q_h*faraday*phi
//     );

//     Info<< "Reading field comp_e\n" << endl;
// volScalarField comp_e
// (
//     IOobject
//     (
//         "comp_e",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     //mesh
//     //interface*exp((mu_e - E_e_P3HT)/(R*T)) + (1.0-interface)*exp((mu_e - E_e_PCBM)/(R*T))
//     T0*exp((mu_e - E_e_P3HT)/(R*T)) + T1*exp((mu_e - E_e_PCBM)/(R*T)) + T2*exp((mu_e - E_e_ter)/(R*T))
// );
// 
//     Info<< "Reading field comp_h\n" << endl;
// volScalarField comp_h
// (
//     IOobject
//     (
//         "comp_h",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     //mesh
//     //interface*c_h_P3HT +  (1.0 - interface)*c_h_PCBM
//     //interface*exp((mu_h - E_h_P3HT)/(R*T)) + (1.0-interface)*exp((mu_h - E_h_PCBM)/(R*T))
//     T0*exp((mu_h - E_h_P3HT)/(R*T)) + T1*exp((mu_h - E_h_PCBM)/(R*T)) + T2*exp((mu_h - E_h_ter)/(R*T)) 
// );
// 
// Info<< "Reading field background\n" << endl;
// volScalarField background
// (
//     IOobject
//     (
//         "background",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     //mesh
//     -1.0*(q_e*comp_e + q_h*comp_h)
// );
// 
// Info<< "Reading field rho\n" << endl;
// volScalarField rho
// (
//     IOobject
//     (
//         "rho",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     //mesh
//     (q_e*comp_e + q_h*comp_h + background)*faraday
// );

//======================= BEGIN ADDED FROM BUNTY =================////
/* CONFLICT: sig1 is declared twice
 * 
 * 
 */

Info<< "Reading smoothening constants\n" << endl;
IOdictionary constProperties
    (
        IOobject
        (
            "constProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar gamma_0
    (
        constProperties.lookup("gamma_0")
    );    
    dimensionedScalar gamma_1
    (
        constProperties.lookup("gamma_1")
    );    
    dimensionedScalar gamma_2
    (
        constProperties.lookup("gamma_2")
    );    
    dimensionedScalar gamma_3
    (
        constProperties.lookup("gamma_3")
    );        
    dimensionedScalar sigt
    (
	constProperties.lookup("sigt")
    );
    dimensionedScalar Tau
    (
        constProperties.lookup("Tau")
    );
    dimensionedScalar Epsilon
    (
        constProperties.lookup("Epsilon")
    );
	
	dimensionedScalar sig
	(
		constProperties.lookup("sig")
	);

	dimensionedScalar sig1
	(
		constProperties.lookup("sig1")
	);

    
//=====================for three variants============================//
    dimensionedScalar DTs
    (
        constProperties.lookup("DTs")
    );
	
	volScalarField lambda_phase
        (
            IOobject
            (
                "lambda_phase",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sigt*((fvc::ddt(T0) + fvc::ddt(T1) + fvc::ddt(T2))/3.0)
        );
	
	volScalarField deltaT0
        (
            IOobject
            (
                "deltaT0",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (sigt*fvc::ddt(T0) - lambda_phase)
        );
	
volScalarField deltaT1
        (
            IOobject
            (
                "deltaT1",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (sigt*fvc::ddt(T1) - lambda_phase)
        );	
	
volScalarField deltaT2
        (
            IOobject
            (
                "deltaT2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (sigt*fvc::ddt(T2) - lambda_phase)
        );	

//======================= END ADDED FROM BUNTY =================////

	