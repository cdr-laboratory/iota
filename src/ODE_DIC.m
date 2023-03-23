function dDICdx = ODE_DIC(x,DIC)

global Rate_org_interp z Kz_water Adv_TH1 R1_carb
  
  Rcarbon = interp1(z,Rate_org_interp,x);  %molCorg/m3/yr
  R1_carb1 = interp1(z,R1_carb,x);  %umol/L/yr
  Adv11_TH = interp1(z,Adv_TH1,x); 
  Kz_water1 = interp1(z,Kz_water,x);
  Oxy1 = -(Rcarbon)*1E3 + R1_carb1 - Adv11_TH; %umol/L/yr
  dDICdx = [ DIC(2) / Kz_water1
           Oxy1];
end