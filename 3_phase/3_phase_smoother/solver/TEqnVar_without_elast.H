//=================== BEGIN Modified double-well potential =======================//
//     volScalarField dwdphi_0  = 18.0*gamma_0*T0*(1.0 -T0)*(1.0 -2.0*T0);
//     volScalarField dwdphi_1  = 18.0*gamma_1*T1*(1.0 -T1)*(1.0 -2.0*T1);
//     volScalarField dwdphi_2  = 18.0*gamma_2*T2*(1.0 -T2)*(1.0 -2.0*T2);
//     volScalarField dwdphi_3  = 18.0*gamma_3*T3*(1.0 -T3)*(1.0 -2.0*T3);
//=================== END Modified double-well potential =========================//

//====================================== BEGIN calculation elastic driving force ==================================================//    
// 	scalar eigen_xx = 0.0; 
// 	scalar eigen_yy = 0.0;
// 	scalar eigen_zz = 0.0;
//     
//         forAll(strain,cellI)
//           {    
// 	    eigen_xx  =   T0[cellI]*T0[cellI]*(3-2*T0[cellI])*cEigenStrain0.component(symmTensor::XX).value();
// 	    eigen_xx +=   T1[cellI]*T1[cellI]*(3-2*T1[cellI])*cEigenStrain1.component(symmTensor::XX).value();
// 	    eigen_xx +=   T2[cellI]*T2[cellI]*(3-2*T2[cellI])*cEigenStrain2.component(symmTensor::XX).value();
// 	    
// 	    eigen_yy  =   T0[cellI]*T0[cellI]*(3-2*T0[cellI])*cEigenStrain0.component(symmTensor::YY).value();
// 	    eigen_yy +=   T1[cellI]*T1[cellI]*(3-2*T1[cellI])*cEigenStrain1.component(symmTensor::YY).value();
// 	    eigen_yy +=   T2[cellI]*T2[cellI]*(3-2*T2[cellI])*cEigenStrain2.component(symmTensor::YY).value();
// 	    
// 	    eigen_zz  =   T0[cellI]*T0[cellI]*(3-2*T0[cellI])*cEigenStrain0.component(symmTensor::ZZ).value();
// 	    eigen_zz +=   T1[cellI]*T1[cellI]*(3-2*T1[cellI])*cEigenStrain1.component(symmTensor::ZZ).value();
// 	    eigen_zz +=   T2[cellI]*T2[cellI]*(3-2*T2[cellI])*cEigenStrain2.component(symmTensor::ZZ).value();
// 	    
//             strain[cellI].component(symmTensor::XX) = gradD[cellI].component(tensor::XX) - eigen_xx;
//             strain[cellI].component(symmTensor::YY) = gradD[cellI].component(tensor::YY) - eigen_yy;
//             strain[cellI].component(symmTensor::ZZ) = gradD[cellI].component(tensor::ZZ) - eigen_zz;
//             strain[cellI].component(symmTensor::XY) = 0;
//             strain[cellI].component(symmTensor::XZ) = 0;
//             strain[cellI].component(symmTensor::YZ) = 0;
//            }
//            
// //    vector eStrain = (T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2);
//            
//    Sigma = (2*(mu2*T3*T3*(3-2*T3) + mu1*(1-T3)*(1-T3)*(1+2*T3))*(symm(fvc::grad(D)) 
// 	      - (T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2)) 
// 	  + (lambda2*T3*T3*(3-2*T3) + lambda1*(1-T3)*(1-T3)*(1+2*T3))*(I*tr(fvc::grad(D) 
// 	      - (T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2))))
// 	  + (mu2_*T3*T3*(3-2*T3) + mu1_*(1-T3)*(1-T3)*(1+2*T3))*strain;		
// 		
//    deltaSigmaD = ((mu1)*twoSymm(fvc::grad(D)) + (lambda1)*(I*tr(fvc::grad(D))) 
// 		- 2*(mu1)*(T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2) 
// 		- (lambda1)*(I*tr((T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2))))
// 		+ (mu1_)*strain;

//====================================== END calculation elastic driving force ==================================================//    


//=============================================== BEGIN variant evolve equations ===============================================//
    
//For T0
    
    fvScalarMatrix T0Eqn
    (
        sigt*Tau*Epsilon*fvm::ddt(T0) -2.0*gamma_0*Epsilon*sig*fvm::laplacian(DTs, T0) == -(18.0*gamma_0*T0*(1.0 -T0)*(1.0 -2.0*T0))/(Epsilon)
    );
//     T0Eqn -= 6.0*T0*(1-T0)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-T0*T0*(3-2*T0)*cEigenStrain0))-(Sigma && cEigenStrain0));

    
//For T1        
    
    fvScalarMatrix T1Eqn
    (
        sigt*Tau*Epsilon*fvm::ddt(T1) -2.0*gamma_1*Epsilon*sig*fvm::laplacian(DTs, T1) == -(18.0*gamma_1*T1*(1.0 -T1)*(1.0 -2.0*T1))/(Epsilon)
    );
//     T1Eqn -= 6.0*T1*(1-T1)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-T1*T1*(3-2*T1)*cEigenStrain1))-(Sigma && cEigenStrain1));

    
//For T2  
    
    fvScalarMatrix T2Eqn
    (
        sigt*Tau*Epsilon*fvm::ddt(T2) -2.0*gamma_2*Epsilon*sig*fvm::laplacian(DTs, T2) == -(18.0*gamma_2*T2*(1.0 -T2)*(1.0 -2.0*T2))/(Epsilon)
    );
//     T2Eqn -= 6.0*T2*(1-T2)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-T2*T2*(3-2*T2)*cEigenStrain2))-(Sigma && cEigenStrain2));    

    
//For T3     
    
//      deltaSigmaD = ((mu2)*twoSymm(fvc::grad(D)) + (lambda2)*(I*tr(fvc::grad(D))) 
// 		   - 2*(mu2)*(T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2) 
// 		   - (lambda2)*(I*tr((T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2))))
// 		   + (mu2_)*strain;

	/*  ============== NOT NEEDED SMOOTHEING =================
    fvScalarMatrix T3Eqn
    (
        sigt*Tau*Epsilon*fvm::ddt(T3) -2.0*gamma_3*Epsilon*sig1*fvm::laplacian(DTs, T3) == -(18.0*gamma_3*T3*(1.0 -T3)*(1.0 -2.0*T3))/(Epsilon)
    );
    ==================== NOT NEEDED SMOOTHEING */
    
//     T3Eqn -= 6.0*T3*(1-T3)*(0.5*(deltaSigmaD && (symm(fvc::grad(D)))));   
    
    T0Eqn.solve();
    T1Eqn.solve();
    T2Eqn.solve();
//     T3Eqn.solve();
 
//     scalar lambda_phase = 0.0;
    scalar numphases = 3.0;   
    
//     lambda_phase = ((T0-T0Old)+(T3-T3Old))/numphases;
//     lambda_phase = sigt*((fvc::ddt(T0) + fvc::ddt(T3))/numphases);
    lambda_phase = sigt*( (fvc::ddt(T0) + fvc::ddt(T1) + fvc::ddt(T2) )/numphases);
    Info<< "Min/max lambda_phase:" << min(lambda_phase()).value() << '\t'<< max(lambda_phase()).value() << endl;
//     lambda_phase = ((fvc::ddt(T0) + fvc::ddt(T1) + fvc::ddt(T2) + fvc::ddt(T3))/numphases);
//     lambda_phase /= numphases;

//=============== Substracting lambda_phase from evolution equations ==============================//    
    deltaT0 = sigt*fvc::ddt(T0) - lambda_phase;
    deltaT1 = sigt*fvc::ddt(T1) - lambda_phase;
    deltaT2 = sigt*fvc::ddt(T2) - lambda_phase;
//     deltaT3 = sigt*fvc::ddt(T3) - lambda_phase;
//     deltaT0 = T0 -T0Old - lambda_phase;
//     deltaT3 = T3 -T3Old - lambda_phase;
//=============== Substracting lambda_phase from evolution equations ==============================//  
    
  //Update the phase-values such that sum upto 1
    T0 -=  lambda_phase*runTime.deltaTValue();
    T1 -=  lambda_phase*runTime.deltaTValue();
    T2 -=  lambda_phase*runTime.deltaTValue();
//     T3 -=  lambda_phase*runTime.deltaTValue();
  ///////////////////////////////////////////////
//     Info<< "Min/max deltaT0:" << min(deltaT0()).value() << '\t'<< max(deltaT0()).value() << endl;
    Info<< "Min/max T0:" << min(T0()).value() << '\t'<< max(T0()).value() << endl;
	Info<< "Min/max T1:" << min(T0()).value() << '\t'<< max(T0()).value() << endl;
	Info<< "Min/max T2:" << min(T0()).value() << '\t'<< max(T0()).value() << endl;
	
    
/* ****************** NOT NEEDED FOR SMOOTHENING ***********************
//==================== BEGIN calculation for volume preservation ==========================//
    scalar lambdaVol_0 = 0.0;
    scalar lambdaVol_1 = 0.0;
    scalar lambdaVol_2 = 0.0;
    
    scalar sumHprime_0 = 0.0;
    scalar sumHprime_1 = 0.0;
    scalar sumHprime_2 = 0.0;
    
    scalar sumRHS0 = 0.0;
    scalar sumRHS1 = 0.0;
    scalar sumRHS2 = 0.0;

//     scalar sumT = 0.0;
//     scalar sum2T = 0.0;
//     scalar sum3T = 0.0;
	
//     forAll(T0, cellI){
// sumRHS0 = gSum(deltaT0());    
    sumRHS0 = gSum(deltaT0());
    sumRHS1 = gSum(deltaT1());
    sumRHS2 = gSum(deltaT2());
    
    
    volScalarField Hphi_0 =  6.0*T0*(1-T0);
    volScalarField Hphi_1 =  6.0*T1*(1-T1);
    volScalarField Hphi_2 =  6.0*T2*(1-T2);
    
    sumHprime_0 = gSum(Hphi_0());
    sumHprime_1 = gSum(Hphi_1());
    sumHprime_2 = gSum(Hphi_2());
    
    
//     forAll(deltaT0, cellI){
//       sumRHS0 += deltaT0[cellI];
//       sumRHS1 += deltaT1[cellI];
//       sumRHS2 += deltaT2[cellI];
//       sumHprime_0 += 6.0*T0[cellI]*(1-T0[cellI]);
//       sumHprime_1 += 6.0*T1[cellI]*(1-T1[cellI]);
//       sumHprime_2 += 6.0*T2[cellI]*(1-T2[cellI]);
// //       sumT += T0[cellI];
// //       sum2T+= pow(T0[cellI],2);
// //       sum3T += T0[cellI]-T0Old[cellI];
//     }
    
//     Info<< " sumRHS0:" << sumRHS0 << endl;
//     Info<< " sumRHS1:" << sumRHS1 << endl;
//     Info<< " sumRHS2:" << sumRHS2 << endl;
    
//      scalar lambdaVol_0 = 2.0*sum3T/(6.0*(sumT-sum2T));
//     lambdaVol_0  = (2.0*sumRHS0)/sumHprime_0;
    lambdaVol_0  = (2.0*sumRHS0 + sumRHS1 + sumRHS2)/sumHprime_0;
//     lambdaVol_0 /= sumHprime_0;
    
    lambdaVol_1  = (sumRHS0 + 2.0*sumRHS1 + sumRHS2)/sumHprime_1;
//     lambdaVol_1 /= sumHprime_1;
    
    lambdaVol_2  = (sumRHS0 + sumRHS1 + 2.0*sumRHS2)/sumHprime_2;
//     lambdaVol_2 /= sumHprime_2;
//     Info<< " lambdaVol_0:" << lambdaVol_0 << endl;
//     Info<< " lambdaVol_1:" << lambdaVol_1 << endl;
//     Info<< " lambdaVol_2:" << lambdaVol_2 << endl;
    
//     volScalarField lambdaVol_3 = ((6.0*T0*(1-T0)*lambdaVol_0)/numphases);
    volScalarField lambdaVol_3 = ((6.0*T0*(1-T0)*lambdaVol_0 + 6.0*T1*(1-T1)*lambdaVol_1 + 6.0*T2*(1-T2)*lambdaVol_2)/numphases);
//     scalar lambdaVol_3 = ((sumHprime_0*lambdaVol_0 + sumHprime_1*lambdaVol_1 + sumHprime_2*lambdaVol_2)/numphases);
//     lambdaVol_3 /= numphases;
//     forAll(T0, cellI){
//       lambdaVol_3[cellI] = ((6.0*T0[cellI]*(1-T0[cellI])*lambdaVol_0 + 6.0*T1[cellI]*(1-T1[cellI])*lambdaVol_1 
// 			      + 6.0*T2[cellI]*(1-T2[cellI])*lambdaVol_2)/numphases);
//     }
    Info<< "Min/max sum_lambda:" << min((6.0*T0*(1-T0)*lambdaVol_0 + 6.0*T1*(1-T1)*lambdaVol_1 + 6.0*T2*(1-T2)*lambdaVol_2)-numphases*lambdaVol_3).value() << '\t'<< max((6.0*T0*(1-T0)*lambdaVol_0 + 6.0*T1*(1-T1)*lambdaVol_1 + 6.0*T2*(1-T2)*lambdaVol_2)-numphases*lambdaVol_3).value() << endl;
//============================ END calculation for volume preservation ==================================//

//========== Substracting lambda terms from evolution equation for volume preserve ==========//
      deltaT0 = -(6.0*T0*(1-T0)*lambdaVol_0 -lambdaVol_3);
      deltaT1 = -(6.0*T1*(1-T1)*lambdaVol_1 -lambdaVol_3);
      deltaT2 = -(6.0*T2*(1-T2)*lambdaVol_2 -lambdaVol_3);
      deltaT3 = lambdaVol_3;

//     forAll(deltaT0, cellI){
//       deltaT0[cellI] = -(6.0*T0[cellI]*(1-T0[cellI])*lambdaVol_0 -lambdaVol_3);
//       deltaT1[cellI] = -(6.0*T1[cellI]*(1-T1[cellI])*lambdaVol_1 -lambdaVol_3);
//       deltaT2[cellI] = -(6.0*T2[cellI]*(1-T2[cellI])*lambdaVol_2 -lambdaVol_3);
//       deltaT3[cellI] = -lambdaVol_3;
//     }
  
//========== Substracting lambda terms from evolution equation for volume preserve ============//
    
//=========================== phase field update =============================//
//     T0 += (deltaT0 - sigt*fvc::ddt(T0))*runTime.deltaTValue();	//precipitate1
//     T1 += (deltaT1 - sigt*fvc::ddt(T1))*runTime.deltaTValue();	//precipitate2
//     T2 += (deltaT2 - sigt*fvc::ddt(T2))*runTime.deltaTValue();	//precipitate3
//     T3 += (deltaT3 - sigt*fvc::ddt(T3))*runTime.deltaTValue();	//matrix

//     forAll(T0, cellI){
//       T0[cellI] += - (6.0*T0[cellI]*(1-T0[cellI])*lambdaVol_0 -lambdaVol_3[cellI]);
//       T3[cellI] += - (lambdaVol_3[cellI]);
//     }
    
    T0 +=  deltaT0*runTime.deltaTValue();
    T1 +=  deltaT1*runTime.deltaTValue();
    T2 +=  deltaT2*runTime.deltaTValue();
    T3 +=  deltaT3*runTime.deltaTValue();
//     forAll(T0, cellI){
// 	T0[cellI] +=  -(6.0*T0[cellI]*(1-T0[cellI])*lambdaVol_0 -lambdaVol_3[cellI]);
// 	T3[cellI] +=  -lambdaVol_3[cellI];
//     }
//=========================== phase field update =============================//
*/

//====================================== END variant evolve equations ==================================================//