function res = BC_DIC(DICa,DICb)

global DICinit

  res = [ DICa(1)-DICinit 
          DICb(2) ];
 
end