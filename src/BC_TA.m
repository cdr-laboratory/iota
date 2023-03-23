function res = BC_TA(TAa,TAb)

global TAinit

  res = [ TAa(1)-TAinit 
          TAb(2) ];
 
end