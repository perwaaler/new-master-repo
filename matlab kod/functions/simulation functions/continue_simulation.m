function out = continue_simulation(A0,B0,xinit,disable_crash)
% used to determine when to end a simulated encounter

if disable_crash==1
  engage_EA = double(imag(A0)<4 || real(B0)>-xinit); 
else
  engage_EA = double(real(A0) < xinit   && imag(A0)<4 && real(B0) > -xinit && imag(A0)<imag(B0)+1.5 && real(A0)<real(B0)+1.5);
end
out = engage_EA;
end