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

function vehicle = aero_analysis(mission, vehicle)

% Add missing segments vector to vehicle components
for i = 1 : length(vehicle.components)
    if (is_type(vehicle.components{i}, 'fuselage') || is_type(vehicle.components{i}, 'wing') || is_type(vehicle.components{i}, 'boom'))
        vehicle.components{i}.segments = repmat({struct()}, length(mission.segments), 1);
        for j = 1 : length(vehicle.components{i}.segments)
            vehicle.components{i}.segments{j}.name = mission.segments{j}.name;
            vehicle.components{i}.segments{j}.base_drag_coefficient = 0;
            vehicle.components{i}.segments{j}.lift_slope_coefficient = 0;
        end
    end
end

% Add missing segments vector to vehicle
vehicle.segments = repmat({struct()}, length(mission.segments), 1);
for i = 1 : length(vehicle.segments)
    vehicle.segments{i}.name = mission.segments{i}.name;
    vehicle.segments{i}.base_drag_coefficient = 0;
end

% Base drag coefficient calculations
for i = 1 : length(vehicle.components)
    if is_type(vehicle.components{i}, 'fuselage')
        for j = 1 : length(mission.segments)
            [~, comp_segment_id] = find_by_name(vehicle.components{i}.segments, mission.segments{j}.name);

            vehicle.components{i}.segments{comp_segment_id}.base_drag_coefficient = friction_coeff(vehicle.components{i}.length, mean(abs(mission.segments{j}.velocity)), mean(mission.segments{j}.speed_sound), mean(mission.segments{j}.density), air_viscosity(mean(mission.segments{j}.temperature))) *...
                fuselage_form_factor(vehicle.components{i}.length, vehicle.components{i}.diameter) *...
                vehicle.components{i}.interf_factor *...
                vehicle.components{i}.area_wet / vehicle.area_ref;

            [~, vehicle_segment_id] = find_by_name(vehicle.segments, mission.segments{j}.name);

            vehicle.segments{vehicle_segment_id}.base_drag_coefficient = vehicle.segments{vehicle_segment_id}.base_drag_coefficient + vehicle.components{i}.segments{comp_segment_id}.base_drag_coefficient;
        end
    elseif is_type(vehicle.components{i}, 'wing')
        for j = 1 : length(mission.segments)
            [~, comp_segment_id] = find_by_name(vehicle.components{i}.segments, mission.segments{j}.name);

            m = mean(abs(mission.segments{j}.velocity)) / mean(mission.segments{j}.speed_sound);
            bb = sqrt(1 - m^2);
            
            if is_type(vehicle.components{i}, 'wing.vtail')
                vehicle.components{i}.segments{comp_segment_id}.lift_slope_coefficient = vehicle.components{i}.airfoil.lift_slope_coefficient * vehicle.components{i}.aspect_ratio /...
                (2 + sqrt(4 + vehicle.components{i}.aspect_ratio^2 * bb^2 * (1 + tand(sweep_xcmax(vehicle.components{i},vehicle.components{i}.airfoil.xc_max,vehicle.components{i}.root_chord,vehicle.components{i}.taper_ratio,vehicle.components{i}.half_span,vehicle.components{i}.sweep_le))^2 / bb^2)));

                vehicle.components{i}.segments{comp_segment_id}.base_drag_coefficient = friction_coeff(vehicle.components{i}.mean_chord, mean(abs(mission.segments{j}.velocity)), mean(mission.segments{j}.speed_sound), mean(mission.segments{j}.density), air_viscosity(mean(mission.segments{j}.temperature))) *...
                wing_form_factor(vehicle.components{i}.airfoil.xc_max, vehicle.components{i}.airfoil.tc_max, sweep_xcmax(vehicle.components{i},vehicle.components{i}.airfoil.xc_max,vehicle.components{i}.root_chord,vehicle.components{i}.taper_ratio,vehicle.components{i}.half_span,vehicle.components{i}.sweep_le), m) *...
                vehicle.components{i}.interf_factor *...
                vehicle.components{i}.area_wet / vehicle.area_ref;
            else
                vehicle.components{i}.segments{comp_segment_id}.lift_slope_coefficient = vehicle.components{i}.airfoil.lift_slope_coefficient * vehicle.components{i}.aspect_ratio /...
                (2 + sqrt(4 + vehicle.components{i}.aspect_ratio^2 * bb^2 * (1 + tand(sweep_xcmax(vehicle.components{i},vehicle.components{i}.airfoil.xc_max,vehicle.components{i}.root_chord,vehicle.components{i}.taper_ratio,vehicle.components{i}.span,vehicle.components{i}.sweep_le))^2 / bb^2)));

                vehicle.components{i}.segments{comp_segment_id}.base_drag_coefficient = friction_coeff(vehicle.components{i}.mean_chord, mean(abs(mission.segments{j}.velocity)), mean(mission.segments{j}.speed_sound), mean(mission.segments{j}.density), air_viscosity(mean(mission.segments{j}.temperature))) *...
                wing_form_factor(vehicle.components{i}.airfoil.xc_max, vehicle.components{i}.airfoil.tc_max, sweep_xcmax(vehicle.components{i},vehicle.components{i}.airfoil.xc_max,vehicle.components{i}.root_chord,vehicle.components{i}.taper_ratio,vehicle.components{i}.span,vehicle.components{i}.sweep_le), m) *...
                vehicle.components{i}.interf_factor *...
                vehicle.components{i}.area_wet / vehicle.area_ref;
            end
            
            [~, vehicle_segment_id] = find_by_name(vehicle.segments, mission.segments{j}.name);

            vehicle.segments{vehicle_segment_id}.base_drag_coefficient = vehicle.segments{vehicle_segment_id}.base_drag_coefficient + vehicle.components{i}.segments{comp_segment_id}.base_drag_coefficient;
        end
    elseif is_type(vehicle.components{i}, 'boom')
        for j = 1 : length(mission.segments)
            [~, comp_segment_id] = find_by_name(vehicle.components{i}.segments, mission.segments{j}.name);
            [~, vehicle_segment_id] = find_by_name(vehicle.segments, mission.segments{j}.name);
            for k = 1 : size(vehicle.components{i}.position.structure,1)     
                cf = friction_coeff(vehicle.components{i}.geometry.chord(k), mean(abs(mission.segments{j}.velocity))*cosd(vehicle.components{i}.geometry.angle(k)), mean(mission.segments{j}.speed_sound), mean(mission.segments{j}.density), air_viscosity(mean(mission.segments{j}.temperature)));
        
                vehicle.components{i}.segments{comp_segment_id}.base_drag_coefficient = (2*cf*(1+vehicle.components{i}.geometry.tc_max(k))+vehicle.components{i}.geometry.tc_max(k)^2)*vehicle.components{i}.geometry.length(k)*vehicle.components{i}.geometry.chord(k)/vehicle.area_ref; 
        
                vehicle.segments{vehicle_segment_id}.base_drag_coefficient = vehicle.segments{vehicle_segment_id}.base_drag_coefficient + vehicle.components{i}.segments{comp_segment_id}.base_drag_coefficient;
            end
        end
    end
end

function f = sweep_xcmax(wing,xc,root_chord,taper_ratio,span,sweep_le)
if is_type(wing, 'wing.vtail')
    f = atand((xc*root_chord*(taper_ratio-1))/span+tand(sweep_le));
else
    f = atand((2*xc*root_chord*(taper_ratio-1))/span+tand(sweep_le));
end

function f = fuselage_form_factor(l, d)
ld = l / d;
f = 1 + 60 / ld^3 + ld / 400;

function f = wing_form_factor(xc_max, tc_max, sweep_tc_max, m)
f = (1 + 0.6 / xc_max * tc_max + 100 * tc_max^4) * 1.34 * m^0.18 * cosd(sweep_tc_max)^0.28;

