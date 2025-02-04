% Aircraft design tool
%
% Copyright (C) 2022 Mario Bras
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 3 as
% published by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

function vehicle = build_vehicle(vehicle)
vehicle.area_ref = 0;
% Add missing area_ref to fuselage 
for i = 1 : length(vehicle.components)
    if is_type(vehicle.components{i}, 'fuselage')
        vehicle.components{i}.area_wet = fuselage_area_wet(vehicle.components{i}.length, vehicle.components{i}.diameter);
        [~,f] = find_by_type(vehicle.components, 'fuselage');
    elseif is_type(vehicle.components{i}, 'wing')
        vehicle.components{i}.sweep_c4 = sweep_c4(vehicle.components{i}.aspect_ratio,vehicle.components{i}.taper_ratio,vehicle.components{i}.sweep_le);
        if is_type(vehicle.components{i}, 'wing.vtail')
            vehicle.components{i}.root_chord = root_chord(vehicle.components{i}.taper_ratio,2*vehicle.components{i}.half_span, vehicle.components{i}.aspect_ratio);
            vehicle.components{i}.mean_chord = mean_chord(vehicle.components{i}.taper_ratio,vehicle.components{i}.root_chord);
            vehicle.components{i}.area_ref = ((2 * vehicle.components{i}.half_span)^2 / vehicle.components{i}.aspect_ratio)/2;
            vehicle.components{i}.area_wet = wing_area_wet(vehicle.components{i}.airfoil.tc_max, vehicle.components{i}.area_ref);
        else
            vehicle.components{i}.root_chord = root_chord(vehicle.components{i}.taper_ratio,vehicle.components{i}.span, vehicle.components{i}.aspect_ratio);
            vehicle.components{i}.mean_chord = mean_chord(vehicle.components{i}.taper_ratio,vehicle.components{i}.root_chord);
            vehicle.components{i}.area_ref = vehicle.components{i}.span^2 / vehicle.components{i}.aspect_ratio;
            vehicle.components{i}.area_wet = wing_area_wet(vehicle.components{i}.airfoil.tc_max, vehicle.components{i}.area_ref);
        end
        vehicle.components{i}.tip_chord = vehicle.components{i}.root_chord * vehicle.components{i}.taper_ratio;
    elseif is_type(vehicle.components{i}, 'boom') % calculate angle with respect to horizontal line
         for j = 1:size(vehicle.components{i}.position.structure,1)
             delta_x = vehicle.components{i}.position.structure(j,1,1) - vehicle.components{i}.position.structure(j,2,1); % difference of x coordinates between start and end points of boom
             boom_length = sqrt(delta_x^2 + (vehicle.components{i}.position.structure(j,1,2) - vehicle.components{i}.position.structure(j,2,2))^2);
             vehicle.components{i}.geometry.angle(j) = abs(asind(delta_x/boom_length));
         end
    end
    
    % Calculate aircraft's reference area using all lifting surfaces
    if is_type(vehicle.components{i}, 'wing.main') || is_type(vehicle.components{i}, 'wing.canard') || is_type(vehicle.components{i}, 'wing.secondary')
        vehicle.area_ref = vehicle.area_ref + vehicle.components{i}.area_ref;
    end
end

if vehicle.area_ref == 0 % TEMPORARY: using fuselage wetted area as reference for aerodynamic calculations of rotorcraft
    vehicle.area_ref = vehicle.components{f}.area_wet;
end

function a = fuselage_area_wet(l, d)
a = pi() * d * l + pi() * d^2;

function a = wing_area_wet(tc_max, wing_area_ref)
if (tc_max > 0.05)
    a = (1.977 + 0.52 * tc_max) * wing_area_ref;
else
    a = 2.003 * wing_area_ref;
end

function f = root_chord(taper_ratio,span,ar)
f = 2 * span / (ar * (1 + taper_ratio));

function f = mean_chord(taper_ratio,root_chord)
f = 2/3 * root_chord * (1 + taper_ratio + taper_ratio^2) /(1 + taper_ratio);

function f = sweep_c4(AR,taper_ratio,sweep_le)
f = atand(tand(sweep_le)-(1-taper_ratio)/(AR*(1+taper_ratio)));
