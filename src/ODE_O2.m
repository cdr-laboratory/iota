function dpdx = ODE_O2(x,p)

global Rate_org_interp z kO2 Kz_water Iron kFe Adv_TH
  
  Rcarbon = interp1(z,Rate_org_interp,x);  %molCorg/m3/yr
  Adv1_TH = interp1(z,Adv_TH,x);
  Fe_1 = interp1(z,Iron,x);
  Kz_water1 = interp1(z,Kz_water,x);
  Oxy = (Rcarbon.*p(1)/(p(1)+kO2))*1E3 + kFe.*Fe_1.*p(1) - Adv1_TH; %umol/L/yr
  dpdx = [ p(2) / Kz_water1
           Oxy];
end