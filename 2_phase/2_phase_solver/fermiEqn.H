E_e_P3HT = -q_e*faraday*eqm_pot_diff - R*T*log(c_e_P3HT/c_e_PCBM);             // Fermi Energies
E_e_PCBM =  0.0;
E_h_P3HT = -q_h*faraday*eqm_pot_diff - R*T*log(c_h_P3HT/c_h_PCBM);             
E_h_PCBM =  0.0;

Info  << "  E_e_P3HT = " << E_e_P3HT.value() << endl;
Info  << "  E_e_PCBM = " << E_e_PCBM.value() << endl;
Info  << "  E_h_P3HT = " << E_h_P3HT.value() << endl;
Info  << "  E_h_PCBM = " << E_h_PCBM.value() << nl << endl;
    
