function dipolestat1D=DipoleCalculation1D(paux,PSI,x, y, z,dim)
Current=cd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO:
% The function calculate the dynamic dipole for a given atomic species specified
% by the parameter n l m Z
% For all the value of the velocity vector v_au
%
% Assumed 1D interaction z-polarization.
% SFEXISS F.V 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Evaluate the dipole between the state with momentum p=v-A and the
% bounded state
U=length(paux);
dx=dim(2)-dim(1);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dipolestat1D = zeros(1,U);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static Dipole
% Assuming z-oriented 

f = waitbar(0/U,'Calculating Static Dipole...');   
for i=1:U
waitbar(i/U)   
Ip= 1/(2*pi)^(3/2)*PSI.*exp(-1i*paux(i).*z);
%dipolestatA(1,i) = sum(sum(sum( x .* Ip )))*dx.^3;
%dipolestatA(2,i) = sum(sum(sum( y .* Ip )))*dx.^3;
dipolestat1D(i) = sum(sum(sum( z .* Ip )))*dx.^3;
end

disp("Static Dipole Calculated")
disp(" ");
