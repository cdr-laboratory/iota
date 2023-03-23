function dTAdx = ODE_TA(x,TA)

global z Kz_water R1_carb R_diss_mineral_1 Adv_TH2 
  
  R1_carb1 = interp1(z,R1_carb,x);  %umol/L/yr
  R1_diss_mineral_1 = interp1(z,R_diss_mineral_1,x);  %umol/L/yr
  Adv11_TH = interp1(z,Adv_TH2,x); 
  Kz_water1 = interp1(z,Kz_water,x);
  NR = + 2.*R1_carb1 - R1_diss_mineral_1 - Adv11_TH; %umol/L/yr
  dTAdx = [ TA(2) / Kz_water1
           NR];
end