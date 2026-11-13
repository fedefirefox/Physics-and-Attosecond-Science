function dipoleAux = dipoleShift(dipolestat1D,Asim_IR,paux)
% d(p)--->d(p+A_IR) I want to shift the dipole of the quantity p

% plot(paux,abs(dipolestat1D),'ob');
% hold on
% plot(paux,unwrap(angle(dipolestat1D)),'ob');

pauxNew=linspace(paux(1),paux(end),10e4);
CrossSection=abs(dipolestat1D);
Phase=unwrap(angle(dipolestat1D));

CrossSectionNew=interp1(paux,CrossSection,pauxNew,'nearest');
PhaseNEw=interp1(paux,Phase,pauxNew,'spline');

% 
% 
% plot(pauxNew,abs(dipolestat1DNew),'v');
% hold on
% plot(pauxNew,unwrap(angle(dipolestat1DNew)),'v');
% 

shiftAmount=Asim_IR/(abs(pauxNew(2)-pauxNew(1)));

shiftAmount=round(shiftAmount);
%shiftAmount=0;

CrossSectionNew=circshift(CrossSectionNew,shiftAmount);
PhaseNEw=circshift(PhaseNEw,shiftAmount);
if shiftAmount<0
CrossSectionNew(end+shiftAmount:end)=CrossSectionNew(end);
PhaseNEw(end+shiftAmount:end)=PhaseNEw(end);
elseif shiftAmount>0
CrossSectionNew(1:shiftAmount)=CrossSectionNew(1);    
PhaseNEw(1:shiftAmount)=PhaseNEw(1);       
    else
        0;
end
% 
% plot(pauxNew,abs(dipolestat1DNew),'vk');
% plot(pauxNew,unwrap(angle(dipolestat1DNew)),'vk');


CrossAux=smooth(interp1(pauxNew,CrossSectionNew,paux,'nearest'));
PhaseAux=smooth(interp1(pauxNew,PhaseNEw,paux,'spline'));

% 
% plot(paux,CrossAux,'r');
% hold on
% plot(paux,PhaseAux,'r');
% hold off
% pause(0.01)

dipoleAux=CrossAux.*exp(1i*PhaseAux);

end

