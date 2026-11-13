function [signal, amplitude] = sfexis_spectrum(time, xuv, potential, energy, ip)
%SFEXIS_SPECTRUM Constant-dipole SFA integral, all inputs in atomic units.
% time and potential: Nt samples; xuv: Nt x Ndelay; energy: Ne samples.
% Returns Ndelay x Ne arrays. The p weighting preserves the legacy convention
% for a directional spectrum from a 3D momentum density (see docs/GUIDE.md).
% Original SFEXIS: Federico Vismarra. Cite doi:10.1088/2515-7647/ac7991.
validateattributes(time, {'double'}, {'vector','real','finite','nonempty'});
validateattributes(potential, {'double'}, {'vector','real','finite','nonempty'});
validateattributes(xuv, {'double'}, {'2d','finite','nonempty'});
validateattributes(energy, {'double'}, {'vector','real','finite','nonnegative','nonempty'});
validateattributes(ip, {'double'}, {'scalar','real','finite','positive'});
time = time(:);
potential = potential(:);
energy = energy(:)';
if numel(time) < 3 || any(diff(time) <= 0) || ...
        numel(potential) ~= numel(time) || size(xuv,1) ~= numel(time) || ...
        any(diff(energy) <= 0)
    error('sfexis:InvalidGrid', 'Check increasing axes and time dimensions.');
end
amplitude = complex(zeros(size(xuv,2), numel(energy)));
momentum = sqrt(2*energy);
for k = 1:numel(energy)
    % Integral from t to t_end. Reverse the grid, then restore its sign.
    integrand = momentum(k)*potential + potential.^2/2;
    tail = -flipud(cumtrapz(flipud(time), flipud(integrand)));
    phase = exp(1i*((energy(k)+ip)*time - tail));
    amplitude(:,k) = (-1i*trapz(time, xuv.*phase, 1)).';
end
signal = abs(amplitude).^2 .* momentum;
end
