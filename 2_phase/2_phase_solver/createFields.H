Info<< "Reading physicalProperties\n" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar epsilon0
(
    physicalProperties.lookup("epsilon0")
);

dimensionedScalar k
(
    physicalProperties.lookup("k")
);

dimensionedScalar q_e
(
    physicalProperties.lookup("q_e")
);

dimensionedScalar q_h
(
    physicalProperties.lookup("q_h")
);

dimensionedScalar dielectric
(
    physicalProperties.lookup("dielectric")
);

dimensionedScalar permittivity
(
    physicalProperties.lookup("permittivity")
);

dimensionedScalar R
(
    physicalProperties.lookup("R")
);

dimensionedScalar T
(
    physicalProperties.lookup("T")
);

dimensionedScalar faraday
(
    physicalProperties.lookup("faraday")
);

dimensionedScalar D_e_P3HT
(
    physicalProperties.lookup("D_e_P3HT")
);

dimensionedScalar D_e_PCBM
(
    physicalProperties.lookup("D_e_PCBM")
);

    dimensionedScalar D_h_P3HT
(
    physicalProperties.lookup("D_h_P3HT")
);

dimensionedScalar D_h_PCBM
(
    physicalProperties.lookup("D_h_PCBM")
);
dimensionedScalar c_h_P3HT
(
    physicalProperties.lookup("c_h_P3HT")
);

dimensionedScalar c_e_P3HT
(
    physicalProperties.lookup("c_e_P3HT")
);
dimensionedScalar c_h_PCBM
(
    physicalProperties.lookup("c_h_PCBM")
);
dimensionedScalar c_e_PCBM
(
    physicalProperties.lookup("c_e_PCBM")
);
dimensionedScalar G0
(
    physicalProperties.lookup("G0")
);
dimensionedScalar zeta
(
    physicalProperties.lookup("zeta")
);
dimensionedScalar R0
(
    physicalProperties.lookup("R0")
);
dimensionedScalar eqm_pot_diff
(
    physicalProperties.lookup("eqm_pot_diff")
);
dimensionedScalar E_e_P3HT
(
    physicalProperties.lookup("E_e_P3HT")
);
dimensionedScalar E_h_P3HT
(
    physicalProperties.lookup("E_h_P3HT")
);
dimensionedScalar E_e_PCBM
(
    physicalProperties.lookup("E_e_PCBM")
);
dimensionedScalar E_h_PCBM
(
    physicalProperties.lookup("E_h_PCBM")
);
dimensionedScalar sig
(
    physicalProperties.lookup("sig")
);
dimensionedScalar sig1
(
    physicalProperties.lookup("sig1")
);

volScalarField interface
(
    IOobject
    (
        "interface",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field phi\n" << endl;
volScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField mu_e
(
    IOobject
    (
        "mu_e",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //(E_e_P3HT+R*T*log(comp_e))*interface+(1-interface)*E_e_PCBM+R*T*log(comp_e)
);

volScalarField mu_h
(
    IOobject
    (
        "mu_h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //(E_h_P3HT+R*T*log(comp_h))*interface+(1-interface)*E_h_PCBM+R*T*log(comp_h)
);

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

    Info<< "Reading field comp_e\n" << endl;
volScalarField comp_e
(
    IOobject
    (
        "comp_e",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //interface*c_e_P3HT +  (1.0 - interface)*c_e_PCBM
    //mesh
    interface*exp((mu_e - E_e_P3HT)/(R*T)) + (1.0-interface)*exp((mu_e - E_e_PCBM)/(R*T))
);

    Info<< "Reading field comp_h\n" << endl;
volScalarField comp_h
(
    IOobject
    (
        "comp_h",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //mesh
    //interface*c_h_P3HT +  (1.0 - interface)*c_h_PCBM
    interface*exp((mu_h - E_h_P3HT)/(R*T)) + (1.0-interface)*exp((mu_h - E_h_PCBM)/(R*T))
);

Info<< "Reading field background\n" << endl;
volScalarField background
(
    IOobject
    (
        "background",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    //mesh
    -1.0*(q_e*comp_e + q_h*comp_h)
);

Info<< "Reading field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //mesh
    (q_e*comp_e + q_h*comp_h + background)*faraday
);