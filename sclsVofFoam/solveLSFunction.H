// solve Level-Set function as the re-initialization equation
   Info<< "solve the reinitialization equation"     
       << nl << endl;

   psi == psi0;

   for (int corr=0; corr<int(epsilon.value()/deltaTau.value()); corr++)
   {
      psi = psi + psi0/mag(psi0)*(double(1)-mag(fvc::grad(psi)*dimChange))*deltaTau;
      psi.correctBoundaryConditions();
   }


// update Dirac function
   forAll(mesh.cells(),celli)
   {
      if(mag(psi[celli]) > epsilon.value())
         delta[celli] = double(0);
      else
         delta[celli] = double(1.0)/(double(2.0)*epsilon.value())*(double(1.0)+Foam::cos(M_PI*psi[celli]/epsilon.value()));
   };

// update Heaviside function

   forAll(mesh.cells(),celli)
   {
      if(psi[celli] < -epsilon.value())
         H[celli] = double(0);
      else if(epsilon.value() < psi[celli])
         H [celli] = double(1);
      else
         H[celli] = double(1.0)/double(2.0)*(double(1.0)+psi[celli]/epsilon.value()+Foam::sin(M_PI*psi[celli]/epsilon.value())/M_PI);
   };
//alpha1 = H;
//alpha2 = 1.0 - alpha1;

 const volScalarField limitedH
    (
        "limitedH",
        min(max(H, scalar(0)), scalar(1))
    );

rho == rho2 + (rho1 - rho2)*limitedH;
//rho == (rho2 + (rho1 - rho2)*(Foam::exp(25.00*limitedH)-1.00)/Foam::exp(25.00));
//rho == (rho2 + (rho1 - rho2)*pow(limitedH, 1.5*(nuMod1/nuMod2) + 0.75));

volScalarField& nuTemp = const_cast<volScalarField&>(mixture.nu()());
nuTemp == (rho2*nuMod2 + (rho1*nuMod1 - rho2*nuMod2)*limitedH)/(rho2 + (rho1 - rho2)*limitedH);
//nuTemp == (rho2*nuMod2 + (rho1*nuMod1 - rho2*nuMod2)*((Foam::exp(25.00*limitedH)-1.00)/Foam::exp(25.00)))/(rho2 + (rho1 - rho2)*(Foam::exp(25.00*limitedH)-1.00)/Foam::exp(25.00));
//nuTemp == (rho2*nuMod2 + (rho1*nuMod1 - rho2*nuMod2)*pow(limitedH, 1.5*(nuMod1/nuMod2) + 0.75))/(rho2 + (rho1 - rho2)*pow(limitedH, 1.5*(nuMod1/nuMod2) + 0.75));

H == limitedH;