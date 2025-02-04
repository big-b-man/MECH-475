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

function vehicle = design_space_analysis(mission, vehicle, energy)
global constants;

wl = 0:5:2000;
dl = 0:5:10000;
pl = 0:0.0005:0.5;
[plf_grid, wl_grid] = meshgrid(pl, wl);
[plv_grid, dl_grid] = meshgrid(pl, dl);
cf_em = ones(length(wl), length(pl));
cf_ice = ones(length(wl), length(pl));
cf_fc = ones(length(wl), length(pl));
cv = ones(length(dl), length(pl));

% Configure plot
colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'};
figure();
yyaxis left;
cla reset
legend;
hold on;
a = gca;
a.Title.String = 'Design Point';
a.XLim = [0 pl(end)];
a.XLabel.String = 'W/P';
a.YLim = [0 dl(end)];
a.YLabel.String = 'W/A';
a.LineStyleOrder = '-';
colororder(colors)
if is_type(find_by_type(vehicle.components, 'aircraft'),'aircraft.vtol_fixed_wing')
    yyaxis right;
    a.YLim = [0 wl(end)];
    a.YLabel.String = 'W/S';
    a.LineStyleOrder = '-';
    colororder(colors)
end

if is_type(find_by_type(vehicle.components, 'aircraft'),'aircraft.vtol_fixed_wing')
    k = k_parameter(vehicle);

    % Get design wing loading
    wl_design = vehicle.mass * constants.g / vehicle.area_ref;
else
    k = 0;
end

% Get design disk loading and power loading
total_area = 0;
for i = 1 : length(vehicle.components)
    if is_type(vehicle.components{i}, 'driver.prop') && (isfield(vehicle.components{i}.position.cg,'vtol'))
        total_area = total_area + size(vehicle.components{i}.position.cg.vtol,1) * rotor_area(vehicle.components{i});
    end
end
dl_design = vehicle.mass * constants.g / total_area;

fpl_ice_design = 0;
fpl_em_design = 0;
fpl_fc_design = 0;
vpl_design = 0;

% Iterate over flight mission segments
yyaxis left;
forward_region_em = cf_em;
forward_region_ice = cf_ice;
forward_region_fc = cf_fc;
vertical_region = cv;

for i = 1 : length(mission.segments)
    if strcmp(mission.segments{i}.type, 'climb') % Climb segment
        [constraint_ice, constraint_em, vertical_constraint, forward_region_ice, forward_region_em, vertical_region, power_ice, power_em] = climb(plf_grid, plv_grid, wl_grid, dl_grid, wl, dl, k, mission.segments{i}, vehicle, energy);      
        yyaxis left;
        if is_type(find_by_type(vehicle.components, 'aircraft'), 'aircraft.rotary_wing')
            plot(vertical_constraint, dl, 'DisplayName', strcat(mission.segments{i}.name, ": climb constraint - electric motors"));
            vpl_em = vehicle.mass * constants.g / power_em;
            if vpl_em > vpl_design
                vpl_design = vpl_em;
            end
        else
            if isempty(constraint_ice)
                plot(constraint_em, wl, 'DisplayName', strcat(mission.segments{i}.name, ": climb constraint - electric motors"));
                fpl_em = vehicle.mass * constants.g / power_em;
                if fpl_em > fpl_em_design
                    fpl_em_design = fpl_em;
                end
            elseif isempty(constraint_em)
                plot(constraint_ice, wl, 'DisplayName', strcat(mission.segments{i}.name, ": climb constraint - engine"));
                fpl_ice = vehicle.mass * constants.g / power_ice;
                if fpl_ice > fpl_ice_design
                    fpl_ice_design = fpl_ice;
                end
            elseif (~isempty(constraint_ice)) && (~isempty(constraint_em))
                plot(constraint_em, wl, 'DisplayName', strcat(mission.segments{i}.name, ": climb constraint - electric motors"));
                fpl_ice = vehicle.mass * constants.g / power_ice;
                fpl_em = vehicle.mass * constants.g / power_em;
                plot(constraint_ice, wl, 'DisplayName', strcat(mission.segments{i}.name, ": climb constraint - engine"));
                if fpl_ice > fpl_ice_design
                    fpl_ice_design = fpl_ice;
                end
                if fpl_em > fpl_em_design
                    fpl_em_design = fpl_em;
                end
            end
        end
    elseif strcmp(mission.segments{i}.type, 'cruise') % Cruise segment
        [range_constraint, cruise_speed_constraint_ice, cruise_speed_constraint_em , stall_speed_constraint, vertical_constraint, forward_region_ice, forward_region_em, vertical_region, power_ice, power_em] = cruise(plf_grid, plv_grid, wl_grid, dl_grid, wl, dl, k, mission.segments{i}, vehicle, energy);
        yyaxis left;
        
        if is_type(find_by_type(vehicle.components,'aircraft'),'aircraft.rotary_wing')
            plot(vertical_constraint, dl, 'DisplayName', strcat(mission.segments{i}.name, ": forward flight constraint - electric motors"));
                vpl_em = vehicle.mass * constants.g / power_em;
                if vpl_em > vpl_design
                    vpl_design = vpl_em;
                end
        else  
            plot([pl(1) pl(end)], [range_constraint range_constraint], 'DisplayName', strcat(mission.segments{i}.name, ": range constraint - cruise"));
            plot([pl(1) pl(end)], [stall_speed_constraint stall_speed_constraint], 'DisplayName', strcat(mission.segments{i}.name, ": stall speed constraint - cruise"));

            if isempty(cruise_speed_constraint_ice)
                plot(cruise_speed_constraint_em, wl, 'DisplayName', strcat(mission.segments{i}.name, ": cruise speed constraint - electric motors"));
                fpl_em = vehicle.mass * constants.g / power_em;
                if fpl_em > fpl_em_design
                    fpl_em_design = fpl_em;
                end
            elseif isempty(cruise_speed_constraint_em)
                plot(cruise_speed_constraint_ice, wl, 'DisplayName', strcat(mission.segments{i}.name, ": cruise speed constraint - engine"));
                fpl_ice = vehicle.mass * constants.g / power_ice;
                if fpl_ice > fpl_ice_design
                    fpl_ice_design = fpl_ice;
                end
            elseif (~isempty(cruise_speed_constraint_ice)) && (~isempty(cruise_speed_constraint_em))
                plot(cruise_speed_constraint_em, wl, 'DisplayName', strcat(mission.segments{i}.name, ": cruise speed constraint - electric motors"));
                fpl_ice = vehicle.mass * constants.g / power_ice;
                fpl_em = vehicle.mass * constants.g / power_em;
                plot(cruise_speed_constraint_ice, wl, 'DisplayName', strcat(mission.segments{i}.name, ": cruise speed constraint - engine"));
                if fpl_ice > fpl_ice_design
                    fpl_ice_design = fpl_ice;
                end
                if fpl_em > fpl_em_design
                    fpl_em_design = fpl_em;
                end
            end
        end   
    elseif strcmp(mission.segments{i}.type, 'recharge_cruise') % Recharge cruise segment
        [range_constraint, cruise_speed_constraint_ice, cruise_speed_constraint_em , stall_speed_constraint, forward_region_ice, forward_region_em, power_ice, power_em] = recharge_cruise(plf_grid, wl_grid, wl, k, mission.segments{i}, vehicle, energy);
        yyaxis left;
        plot([pl(1) pl(end)], [range_constraint range_constraint], 'DisplayName', strcat(mission.segments{i}.name, ": range constraint - recharge cruise"));
        plot([pl(1) pl(end)], [stall_speed_constraint stall_speed_constraint], 'DisplayName', strcat(mission.segments{i}.name, ": stall speed constraint - recharge cruise"));
        
        if isempty(cruise_speed_constraint_ice)
            plot(cruise_speed_constraint_em, wl, 'DisplayName', strcat(mission.segments{i}.name, ": recharge cruise speed constraint - electric motors"));
            fpl_em = vehicle.mass * constants.g / power_em;
            if fpl_em > fpl_em_design
                fpl_em_design = fpl_em;
            end
        elseif isempty(cruise_speed_constraint_em)
            plot(cruise_speed_constraint_ice, wl, 'DisplayName', strcat(mission.segments{i}.name, ": recharge cruise speed constraint engine"));
            fpl_ice = vehicle.mass * constants.g / power_ice;
            if fpl_ice > fpl_ice_design
                fpl_ice_design = fpl_ice;
            end
        elseif (~isempty(cruise_speed_constraint_ice)) && (~isempty(cruise_speed_constraint_em))
            plot(cruise_speed_constraint_em, wl, 'DisplayName', strcat(mission.segments{i}.name, ": recharge cruise speed constraint - electric motors"));
            fpl_ice = vehicle.mass * constants.g / power_ice;
            fpl_em = vehicle.mass * constants.g / power_em;
            plot(cruise_speed_constraint_ice, wl, 'DisplayName', strcat(mission.segments{i}.name, ": recharge cruise speed constraint - engine"));
            if fpl_ice > fpl_ice_design
                fpl_ice_design = fpl_ice;
            end
            if fpl_em > fpl_em_design
                fpl_em_design = fpl_em;
            end
        end

    elseif strcmp(mission.segments{i}.type, 'hybrid_cruise') % Hybrid cruise segment
        [range_constraint, cruise_speed_constraint_ice, cruise_speed_constraint_em , cruise_speed_constraint_fc, stall_speed_constraint, forward_region_ice, forward_region_em, forward_region_fc, power_ice, power_em, power_fc] = hybrid_cruise(plf_grid, wl_grid, wl, k, mission.segments{i}, vehicle, energy);
        yyaxis left;
        plot([pl(1) pl(end)], [range_constraint range_constraint], 'DisplayName', strcat(mission.segments{i}.name, ": range constraint - hybrid cruise"));
        plot([pl(1) pl(end)], [stall_speed_constraint stall_speed_constraint], 'DisplayName', strcat(mission.segments{i}.name, ": stall speed constraint - hybrid cruise"));
        
        if isempty(cruise_speed_constraint_ice)
            if isempty(cruise_speed_constraint_fc)
                plot(cruise_speed_constraint_em, wl, 'DisplayName', strcat(mission.segments{i}.name, ": hybrid cruise speed constraint - electric motors"));
                fpl_em = vehicle.mass * constants.g / power_em;
                if fpl_em > fpl_em_design
                    fpl_em_design = fpl_em;
                end
            else
                plot(cruise_speed_constraint_em, wl, 'DisplayName', strcat(mission.segments{i}.name, ": hybrid cruise speed constraint - electric motors"));
                fpl_fc = vehicle.mass * constants.g / power_fc;
                fpl_em = vehicle.mass * constants.g / power_em;
                plot(cruise_speed_constraint_fc, wl, 'DisplayName', strcat(mission.segments{i}.name, ": hybrid cruise speed constraint - fuel cell"));
                if fpl_fc > fpl_fc_design
                    fpl_fc_design = fpl_fc;
                end
                if fpl_em > fpl_em_design
                    fpl_em_design = fpl_em;
                end 
            end
        elseif isempty(cruise_speed_constraint_em)
            plot(cruise_speed_constraint_ice, wl, 'DisplayName', strcat(mission.segments{i}.name, ": hybrid cruise speed constraint - engine"));
            fpl_ice = vehicle.mass * constants.g / power_ice;
            if fpl_ice > fpl_ice_design
                fpl_ice_design = fpl_ice;
            end
        elseif (~isempty(cruise_speed_constraint_ice)) && (~isempty(cruise_speed_constraint_em))
            plot(cruise_speed_constraint_em, wl, 'DisplayName', strcat(mission.segments{i}.name, ": hybrid cruise speed constraint - electric motors"));
            fpl_ice = vehicle.mass * constants.g / power_ice;
            fpl_em = vehicle.mass * constants.g / power_em;
            plot(cruise_speed_constraint_ice, wl, 'DisplayName', strcat(mission.segments{i}.name, ": hybrid cruise speed constraint - engine"));
            if fpl_ice > fpl_ice_design
                fpl_ice_design = fpl_ice;
            end
            if fpl_em > fpl_em_design
                fpl_em_design = fpl_em;
            end
        end

    elseif strcmp(mission.segments{i}.type, 'hold') % Hold segment
        [constraint, forward_region_em, forward_region_ice, power_em, power_ice] = loiter(wl_grid, k, mission.segments{i}, vehicle, energy);
        yyaxis left;
        plot([pl(1) pl(end)], [constraint constraint], 'DisplayName', strcat(mission.segments{i}.name, ": endurance constraint"));
        
        if isempty(forward_region_ice)
            fpl_em = vehicle.mass * constants.g / power_em;
            if fpl_em > fpl_em_design
                fpl_em_design = fpl_em;
            end
        elseif isempty(forward_region_em)
            fpl_ice = vehicle.mass * constants.g / power_ice;
            if fpl_ice > fpl_ice_design
                fpl_ice_design = fpl_ice;
            end
        elseif (~isempty(forward_region_ice)) && (~isempty(forward_region_em))
            fpl_ice = vehicle.mass * constants.g / power_ice;
            fpl_em = vehicle.mass * constants.g / power_em;
            if fpl_ice > fpl_ice_design
                fpl_ice_design = fpl_ice;
            end
            if fpl_em > fpl_em_design
                fpl_em_design = fpl_em;
            end
        end
        
   elseif strcmp(mission.segments{i}.type, 'hover') % Hover segment
        [constraint, vertical_region, power] = hover(plv_grid, dl_grid, dl, mission.segments{i}, vehicle, energy);
        yyaxis right;
        plot(constraint, dl, 'DisplayName', strcat(mission.segments{i}.name, ": hover constraint"));

        vpl = vehicle.mass * constants.g / power;
        if vpl > vpl_design
            vpl_design = vpl;
        end
    elseif strcmp(mission.segments{i}.type, 'transition') % Transition segment
        [constraint, vertical_region, power] = transition(plv_grid, dl_grid, wl_design, dl, k, mission.segments{i}, mission.segments{i+1}, vehicle, energy);
        yyaxis right;
        plot(constraint, dl, 'DisplayName', strcat(mission.segments{i}.name, ": transition constraint"));

        vpl = vehicle.mass * constants.g / power;
        if vpl > vpl_design
            vpl_design = vpl;
        end
    elseif strcmp(mission.segments{i}.type, 'vertical_climb') % Vertical climb segment
        [constraint, vertical_region, power] = vertical_climb(plv_grid, dl_grid, dl, mission.segments{i}, vehicle, energy);
        yyaxis right;
        plot(constraint, dl, 'DisplayName', strcat(mission.segments{i}.name, ": vertical climb constraint"));

        vpl = vehicle.mass * constants.g / power;
        if vpl > vpl_design
            vpl_design = vpl;
        end
%     elseif strcmp(mission.segments{i}.type, 'vertical_descent') % Vertical descent segment
%         network = find_network_components(vehicle, find_by_name(energy.networks, mission.segments{i}.energy_network(1)));
% 
%         % TODO
% 
%         vpl = vehicle.mass * constants.g / network_max_power(network);
%         if vpl > vpl_design
%             vpl_design = vpl;
%         end
    end
    
    if (~isempty(forward_region_em))
        cf_em = cf_em .* forward_region_em;
    end
    if (~isempty(forward_region_ice))
        cf_ice = cf_ice .* forward_region_ice;
    end
    if (~isempty(forward_region_fc))
        cf_fc = cf_fc .* forward_region_fc;
    end
    cv = cv .* vertical_region;
end

% Plot feasible design regions
cf_em(~cf_em) = NaN;
cf_ice(~cf_ice) = NaN;
cf_fc(~cf_fc) = NaN;
yyaxis left;
if fpl_em_design ~= 0
    surf(pl, wl, cf_em, 'FaceAlpha', 0.2, 'FaceColor', '#0072BD', 'EdgeColor', 'none', 'DisplayName', 'Forward Flight Design Space Electric Motors');
end
if fpl_ice_design ~= 0
    surf(pl, wl, cf_ice, 'FaceAlpha', 0.2, 'FaceColor', '#77AC30', 'EdgeColor', 'none', 'DisplayName', 'Forward Flight Design Space ICE');
end
if fpl_fc_design ~= 0
    surf(pl, wl, cf_fc, 'FaceAlpha', 0.2, 'FaceColor', '#77AC30', 'EdgeColor', 'none', 'DisplayName', 'Forward Flight Design Space Fuel Cell');
end
cv(~cv) = NaN;
yyaxis right;
surf(pl, dl, cv, 'FaceAlpha', 0.2, 'FaceColor', '#D95319', 'EdgeColor', 'none', 'DisplayName', 'Vertical Flight Design Space');

% Plot design points
yyaxis left;
if fpl_ice_design ~= 0
    scatter(fpl_ice_design, wl_design, 'filled', 'MarkerEdgeColor', '#77AC30', 'MarkerFaceColor', '#77AC30', 'DisplayName', 'Forward Flight Design Point ICE');
end

if fpl_em_design ~= 0
    scatter(fpl_em_design, wl_design, 'filled', 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#0072BD', 'DisplayName', 'Forward Flight Design Point Electric Motors');
end

if fpl_fc_design ~= 0
    scatter(fpl_fc_design, wl_design, 'filled', 'MarkerEdgeColor', '#77AC30', 'MarkerFaceColor', '#77AC30', 'DisplayName', 'Forward Flight Design Point Fuel Cell');
end

yyaxis right;
scatter(vpl_design, dl_design, 'filled', 'MarkerEdgeColor', '#D95319', 'MarkerFaceColor', '#D95319', 'DisplayName', 'Vertical Flight Design Point');

% Helper functions
function [constraint, region, power] = hover(plv_grid, dl_grid, dl, segment, vehicle, energy)
network = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
rotor = find_by_type(network, 'driver.prop');
motor = find_by_type(network, 'motor.prop');

constraint = hover_constraint(dl, segment.density, segment.weight_init, vehicle.mass, rotor{1,1}.tip_velocity, rotor{1,1}.rotor_solidity, rotor{1,1}.base_drag_coefficient, rotor{1,1}.induced_power_factor, motor{1,1}.rating,rotor{1,1}.coaxial,rotor{1,1}.ducted);
region = hover_region(plv_grid, dl_grid, segment.density, segment.weight_init, vehicle.mass, rotor{1,1}.tip_velocity, rotor{1,1}.rotor_solidity, rotor{1,1}.base_drag_coefficient, rotor{1,1}.induced_power_factor, motor{1,1}.rating,rotor{1,1}.coaxial,rotor{1,1}.ducted);

power = network_power(network, 'motor.prop');

% function [constraint, region, power] = transition(plv_grid, dl_grid, wl, dl, k, segment, next_segment, vehicle, energy)
% network = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
% [segment_props, ~] = find_by_name(vehicle.segments, segment.name);
% rotor = find_by_type(network, 'driver.prop');
% 
% constraint = transition_constraint(wl, dl, segment.density, k, segment_props.base_drag_coefficient, rotor.tip_velocity, rotor.rotor_solidity, rotor.base_drag_coefficient, rotor.induced_power_factor, next_segment.velocity, segment.transition_angle);
% region = transition_region(plv_grid, wl, dl_grid, segment.density, k, segment_props.base_drag_coefficient, rotor.tip_velocity, rotor.rotor_solidity, rotor.base_drag_coefficient, rotor.induced_power_factor, next_segment.velocity, segment.transition_angle);
% 
% power = network_power(network, 'motor.prop');

function [constraint, region, power] = vertical_climb(plv_grid, dl_grid, dl, segment, vehicle, energy)
network = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
rotor = find_by_type(network, 'driver.prop');
motor = find_by_type(network, 'motor.prop');

constraint = vertical_climb_constraint(dl, segment.density(1), segment.weight_init, vehicle.mass, rotor{1,1}.tip_velocity, rotor{1,1}.rotor_solidity, rotor{1,1}.base_drag_coefficient, rotor{1,1}.induced_power_factor, segment.velocity, motor{1,1}.rating,rotor{1,1}.coaxial,rotor{1,1}.ducted);
region = vertical_climb_region(plv_grid, dl_grid, segment.density(1), segment.weight_init, vehicle.mass, rotor{1,1}.tip_velocity, rotor{1,1}.rotor_solidity, rotor{1,1}.base_drag_coefficient, rotor{1,1}.induced_power_factor, segment.velocity, motor{1,1}.rating,rotor{1,1}.coaxial,rotor{1,1}.ducted);

power = network_power(network, 'motor.prop');

function [constraint_ice, constraint_em, vertical_constraint, region_ice, region_em, vertical_region, power_ice, power_em] = climb(plf_grid, vpl_grid, wl_grid, dl_grid, wl, dl, k, segment, vehicle, energy)
network = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
source = find_by_type(network, 'energy');
engine = find_by_type(network, 'engine');
motor = find_by_type(network, 'motor');
rotor = find_by_type(network, 'driver.prop');
aircraft = find_by_type(vehicle.components, 'aircraft');
[segment_props, ~] = find_by_name(vehicle.segments, segment.name);

if is_type(source{1,1}, 'energy.fuel')
    if isempty(motor)
        power_em = [];
        region_em = [];
        constraint_em = [];
        if is_type(engine, 'engine.jet')
            constraint_ice = climb_constraint_jet(wl, segment.density(1), segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, segment.angle);
            region_ice = climb_region_jet(plf_grid, wl_grid, segment.density(1), segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, segment.angle);

            power_ice = network_power(network, 'engine.jet');
        elseif is_type(engine, 'engine.prop')
            constraint_ice = climb_constraint_prop(wl, segment.density(1), segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, segment.angle, network_efficiency(network,'engine.prop'), engine{1,1}.rating);
            region_ice = climb_region_prop(plf_grid, wl_grid, segment.density(1), segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, segment.angle, network_efficiency(network,'engine.prop'), engine{1,1}.rating);

            power_ice = network_power(network, 'engine.prop');
        end
    else
        constraint_em = climb_constraint_prop(wl, segment.density(1), segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, segment.angle, network_efficiency(network,'motor.prop'), motor{1,1}.rating);
        region_em = climb_region_prop(plf_grid, wl_grid, segment.density(1), segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, segment.angle, network_efficiency(network,'motor.prop'), motor{1,1}.rating);
        constraint_ice = climb_constraint_prop(wl, segment.density(1), segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, segment.angle, network_efficiency(network,'engine.prop'), engine{1,1}.rating);
        region_ice = climb_region_prop(plf_grid, wl_grid, segment.density(1), segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, segment.angle, network_efficiency(network,'engine.prop'), engine{1,1}.rating);
        power_em = network_power(network, 'motor.prop');
        power_ice = network_power(network, 'engine.prop');
    end
elseif is_type(source{1,1}, 'energy.electric')
    if is_type(aircraft, 'aircraft.rotary_wing')
        vertical_constraint = climb_constraint_rotary_prop(dl, aircraft{1}.rotary_wing, segment, vehicle, rotor, network_efficiency(network,'motor.prop'), motor);
        vertical_region = climb_region_rotary_prop(vpl_grid, dl_grid, aircraft{1}.rotary_wing, segment, vehicle, rotor, network_efficiency(network,'motor.prop'), motor);
        constraint_em = [];
        region_em = [];
    else
        constraint_em = climb_constraint_prop(wl, segment.density(1), segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, segment.angle, network_efficiency(network,'motor.prop'), motor{1,1}.rating);
        region_em = climb_region_prop(plf_grid, wl_grid, segment.density(1), segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, segment.angle, network_efficiency(network,'motor.prop'), motor{1,1}.rating);
    end
    constraint_ice = [];
    region_ice = [];
    power_em = network_power(network, 'motor.prop');
    power_ice = [];
end

function [range_constraint, cruise_speed_constraint_ice, cruise_speed_constraint_em , stall_constraint, vertical_constraint, region_ice, region_em, vertical_region, power_ice, power_em] = cruise(plf_grid, vpl_grid, wl_grid, dl_grid, wl, dl, k, segment, vehicle, energy)
network = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
engine = find_by_type(network, 'engine');
source = find_by_type(network, 'energy');
motor = find_by_type(network, 'motor');
rotor = find_by_type(network, 'driver.prop');
aircraft = find_by_type(vehicle.components, 'aircraft');
[segment_props, ~] = find_by_name(vehicle.segments, segment.name);
main_wing = find_by_type(vehicle.components, 'wing.main');

if ~is_type(aircraft, 'aircraft.rotary_wing')
    stall_constraint = stall_speed_constraint(segment.density, segment.velocity_stall, segment.weight_init, vehicle.mass, main_wing{1,1}.airfoil.cl_max);
    stall_region = stall_speed_region(wl_grid, segment.density, segment.velocity_stall, segment.weight_init, vehicle.mass, main_wing{1,1}.airfoil.cl_max);
else
    stall_constraint = [];
    stall_region = [];
end

if is_type(source{1}, 'energy.fuel')
    if isempty(motor)
        cruise_speed_constraint_em = [];
        power_em = [];
        if is_type(engine, 'engine.jet')
            range_constraint = range_constraint_jet(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            range_region = range_region_jet(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            cruise_speed_constraint_ice = cruise_speed_constraint_jet(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            cruise_speed_region_ice = cruise_speed_region_jet(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            power_ice = network_power(network, 'engine.jet');
        elseif is_type(engine, 'engine.prop')
            range_constraint = range_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            range_region = range_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            cruise_speed_constraint_ice = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);
            cruise_speed_region_ice = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);

            power_ice = network_power(network, 'engine.prop');
        end

        region_ice = range_region.* stall_region.* cruise_speed_region_ice;
        region_em = [];
    else
        range_constraint = range_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
        range_region = range_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

        cruise_speed_constraint_em = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
        cruise_speed_region_em = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
        cruise_speed_constraint_ice = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);
        cruise_speed_region_ice = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);

        power_ice = network_power(network, 'engine.prop');
        power_em = network_power(network, 'motor.prop');

        region_em = range_region.* stall_region.* cruise_speed_region_em;
        region_ice = range_region.* stall_region.* cruise_speed_region_ice;
    end
elseif is_type(source{1,1}, 'energy.electric')
    if is_type(aircraft, 'aircraft.rotary_wing')
        vertical_constraint = fwd_flight_constraint_rotary_wing_prop(dl, aircraft{1}.rotary_wing, segment, vehicle, rotor, network_efficiency(network,'motor.prop'), motor);
        vertical_region = fwd_flight_region_rotary_wing_prop(vpl_grid, dl_grid, aircraft{1}.rotary_wing, segment, vehicle, rotor, network_efficiency(network,'motor.prop'), motor);

        range_constraint = [];

        cruise_speed_constraint_em = [];
        cruise_speed_constraint_ice = [];

        region_em = [];
    else
        range_constraint = range_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
        range_region = range_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

        cruise_speed_constraint_em = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
        cruise_speed_region_em = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
        cruise_speed_constraint_ice = [];
       
        region_em = range_region.* stall_region.* cruise_speed_region_em;
    end
    region_ice = [];
    power_ice = [];

    power_em = network_power(network, 'motor.prop');
end


function [range_constraint, cruise_speed_constraint_ice, cruise_speed_constraint_em , stall_constraint, region_ice, region_em, power_ice, power_em] = recharge_cruise(plf_grid, wl_grid, wl, k, segment, vehicle, energy)
network = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network_prop));
engine = find_by_type(network, 'engine');
source = find_by_type(network, 'energy');
motor = find_by_type(network, 'motor');
[segment_props, ~] = find_by_name(vehicle.segments, segment.name);
main_wing = find_by_type(vehicle.components, 'wing.main');

stall_constraint = stall_speed_constraint(segment.density, segment.velocity_stall, segment.weight_init, vehicle.mass, main_wing{1,1}.airfoil.cl_max);
stall_region = stall_speed_region(wl_grid, segment.density, segment.velocity_stall, segment.weight_init, vehicle.mass, main_wing{1,1}.airfoil.cl_max);

if is_type(source{1}, 'energy.fuel')
    if isempty(motor)
        cruise_speed_constraint_em = [];
        power_em = [];
        if is_type(engine, 'engine.jet')
            range_constraint = range_constraint_jet(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            range_region = range_region_jet(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            cruise_speed_constraint_ice = cruise_speed_constraint_jet(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            cruise_speed_region_ice = cruise_speed_region_jet(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            power_ice = network_power(network, 'engine.jet');
        elseif is_type(engine, 'engine.prop')
            range_constraint = range_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            range_region = range_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            cruise_speed_constraint_ice = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);
            cruise_speed_region_ice = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);

            power_ice = network_power(network, 'engine.prop');
        end

        region_ice = range_region.* stall_region.* cruise_speed_region_ice;
        region_em = [];
    else
        range_constraint = range_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
        range_region = range_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

        cruise_speed_constraint_em = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
        cruise_speed_region_em = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
        cruise_speed_constraint_ice = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);
        cruise_speed_region_ice = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);

        power_ice = network_power(network, 'engine.prop');
        power_em = network_power(network, 'motor.prop');

        region_em = range_region.* stall_region.* cruise_speed_region_em;
        region_ice = range_region.* stall_region.* cruise_speed_region_ice;
    end
elseif is_type(source{1,1}, 'energy.electric')
    range_constraint = range_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
    range_region = range_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

    cruise_speed_constraint_em = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
    cruise_speed_region_em = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
    cruise_speed_constraint_ice = [];
       
    power_ice = [];
    power_em = network_power(network, 'motor.prop');

    region_em = range_region.* stall_region.* cruise_speed_region_em;
    region_ice = [];
end
        
function [range_constraint, cruise_speed_constraint_ice, cruise_speed_constraint_em , cruise_speed_constraint_fc, stall_constraint, region_ice, region_em, region_fc, power_ice, power_em, power_fc] = hybrid_cruise(plf_grid, wl_grid, wl, k, segment, vehicle, energy)
network = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network_fuel));
engine = find_by_type(network, 'engine');
source = find_by_type(network, 'energy');
motor = find_by_type(network, 'motor');
fuel_cell = find_by_type(network, 'fuel_cell');
[segment_props, ~] = find_by_name(vehicle.segments, segment.name);
main_wing = find_by_type(vehicle.components, 'wing.main');

stall_constraint = stall_speed_constraint(segment.density, segment.velocity_stall, segment.weight_init, vehicle.mass, main_wing{1,1}.airfoil.cl_max);
stall_region = stall_speed_region(wl_grid, segment.density, segment.velocity_stall, segment.weight_init, vehicle.mass, main_wing{1,1}.airfoil.cl_max);

if is_type(source{1}, 'energy.fuel')
    power_fc = [];
    cruise_speed_constraint_fc = [];
    region_fc = [];
    if isempty(motor)
        cruise_speed_constraint_em = [];
        power_em = [];
        if is_type(engine, 'engine.jet')
            range_constraint = range_constraint_jet(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            range_region = range_region_jet(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            cruise_speed_constraint_ice = cruise_speed_constraint_jet(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            cruise_speed_region_ice = cruise_speed_region_jet(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            power_ice = network_power(network, 'engine.jet');
        elseif is_type(engine, 'engine.prop')
            range_constraint = range_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            range_region = range_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            cruise_speed_constraint_ice = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);
            cruise_speed_region_ice = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);

            power_ice = network_power(network, 'engine.prop');
        end

        region_ice = range_region.* stall_region.* cruise_speed_region_ice;
        region_em = [];
    else
        range_constraint = range_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
        range_region = range_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

        cruise_speed_constraint_em = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
        cruise_speed_region_em = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
        cruise_speed_constraint_ice = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);
        cruise_speed_region_ice = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'engine.prop'), engine{1,1}.rating);

        power_ice = network_power(network, 'engine.prop');
        power_em = network_power(network, 'motor.prop');

        region_em = range_region.* stall_region.* cruise_speed_region_em;
        region_ice = range_region.* stall_region.* cruise_speed_region_ice;
    end
elseif is_type(source{1}, 'energy.hydrogen')
    power_ice = [];
    cruise_speed_constraint_ice = [];
    region_ice = [];

    range_constraint = range_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
    range_region = range_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

    cruise_speed_constraint_em = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
    cruise_speed_region_em = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'motor.prop'), motor{1,1}.rating);
    cruise_speed_constraint_fc = cruise_speed_constraint_prop(wl, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'fuel_cell'), fuel_cell{1,1}.rating);
    cruise_speed_region_fc = cruise_speed_region_prop(plf_grid, wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k, network_efficiency(network, 'fuel_cell'), fuel_cell{1,1}.rating);

    power_fc = network_power(network, 'fuel_cell');
    power_em = network_power(network, 'motor.prop');

    region_em = range_region.* stall_region.* cruise_speed_region_em;
    region_fc = range_region.* stall_region.* cruise_speed_region_fc;
end

function [constraint, region_em, region_ice, power_em, power_ice] = loiter(wl_grid, k, segment, vehicle, energy)
network = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
source = find_by_type(network, 'energy');
engine = find_by_type(network, 'engine');
motor = find_by_type(network, 'motor');
[segment_props, ~] = find_by_name(vehicle.segments, segment.name);

if is_type(source{1}, 'energy.fuel')
    if isempty(motor)
        region_em = [];
        power_em = [];

        if is_type(engine, 'engine.jet')
            constraint = endurance_constraint_jet(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            region_ice = endurance_region_jet(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            power_ice = network_power(network, 'engine.jet');
        elseif is_type(engine, 'engine.prop')
            constraint = endurance_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
            region_ice = endurance_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

            power_ice = network_power(network, 'engine.prop');
        end
    else 
        constraint = endurance_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
        region_em = endurance_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
        region_ice = [];

        power_em = network_power(network, 'motor.prop');
        power_ice = [];
    end
elseif is_type(source{1}, 'energy.electric')
    if ~isempty(engine)
        constraint = endurance_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
        region_em = endurance_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
        region_ice = endurance_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);

        power_em = network_power(network, 'motor.prop');
        power_ice = network_power(network, 'engine.prop');
    else
        constraint = endurance_constraint_prop(segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
        region_em = endurance_region_prop(wl_grid, segment.density, segment.velocity, segment.weight_init, vehicle.mass, segment_props.base_drag_coefficient, k);
        region_ice = [];

        power_em = network_power(network, 'motor.prop');
        power_ice = [];
    end
end

function k = k_parameter(vehicle)
c = find_by_type(vehicle.components, 'wing.main');
k = 1 / pi / c{1,1}.aspect_ratio / c{1,1}.oswald_efficiency;

%% Performance functions
function v = v_min_thrust(wl, rho, k, cd_0)
v = sqrt(2 * wl / rho * sqrt(k / cd_0));

function v = v_min_power(wl, rho, k, cd_0)
v = sqrt(2 * wl / rho * sqrt(k / 3 / cd_0));

function v = v_best_climb_rate_jet(tl, wl, rho, k, cd_0)
v = sqrt(wl / 3 / rho / cd_0 * (1 / tl + sqrt(1 / tl^2 + 12 * cd_0 * k)));

function v = v_best_climb_rate_prop(wl, rho, k, cd_0)
v = v_min_power(wl, rho, k, cd_0);

function v = v_best_climb_angle_jet(wl, rho, k, cd_0)
v = v_min_thrust(wl, rho, k, cd_0);

function v = v_best_climb_angle_prop(wl, rho, k, cd_0)
v = 0.875 * v_best_climb_rate_prop(wl, rho, k, cd_0); % Raymer pp. 466

function c_l = cl_min_thrust(k, cd_0)
c_l = sqrt(cd_0 / k);

function c_l = cl_min_power(k, cd_0)
c_l = sqrt(3 * cd_0 / k);

function cd = cd_min_thrust(cd_0)
cd = 2 * cd_0;

function cd = cd_min_power(cd_0)
cd = 4 * cd_0;

%% Vertical flight constraint functions
function pl = hover_constraint(dl, rho, m_i, mtom, v_tip, ss, cd, k_i, r, coaxial, ducted)
if strcmp(coaxial,'Yes') % affects induced power (rotor's profile power is affected by disk loading)
    pl = 1./ (1 ./ r .* m_i ./ mtom .* (k_i .* 2.* sqrt(2).* sqrt(2 .* dl .* m_i ./ mtom ./ rho) + rho .* v_tip.^3 .* ss .* cd ./ (dl .* m_i ./ mtom .* 8)));
else
    pl = 1./ (1 ./ r .* m_i ./ mtom .* (k_i .* sqrt(2 .* dl .* m_i ./ mtom ./ rho) + rho .* v_tip.^3 .* ss .* cd ./ (dl .* m_i ./ mtom .* 8)));
end

if strcmp(ducted,'Yes')
    pl = pl .* sqrt(2);
end

function pl = vertical_climb_constraint(dl, rho, m_i, mtom, v_tip, ss, cd, k_i, v_y, r, coaxial, ducted) % affects induced power (rotor's profile power is affected by disk loading)
if strcmp(coaxial, 'Yes')
    pl = 1 ./ (1 ./ r .* m_i ./ mtom .*(v_y - k_i .* 2.* sqrt(2).* v_y ./ 2 + k_i .* 2.* sqrt(2).* sqrt(v_y.^2 + 2 .* dl .* m_i ./ mtom ./ rho) ./ 2 + rho .* v_tip.^3 .* ss .* cd ./ (dl .* m_i ./ mtom .* 8)));
else
    pl = 1 ./ (1 ./ r .* m_i ./ mtom .*(v_y - k_i .* v_y ./ 2 + k_i .* sqrt(v_y.^2 + 2 .* dl .* m_i ./ mtom ./ rho) ./ 2 + rho .* v_tip.^3 .* ss .* cd ./ (dl .* m_i ./ mtom .* 8)));
end

if strcmp(ducted, 'Yes')
    pl = pl .* sqrt(2);
end

function pl = vertical_descent_constraint(dl, rho, v_tip, ss, cd, k_i, v_y, v_i, coaxial, ducted) % affects induced power (rotor's profile power is affected by disk loading)
if strcmp(coaxial,'Yes')
    if v_y / v_i <= -2 % If this condition is met, the vertical climb equation is used for descent, else, an empirical equation is employed
        pl = 1 ./ (v_y - k_i ./ 2 .* sqrt(2) .* 2 .* (v_y + sqrt(v_y.^2 - 2 .* dl ./ rho)) + rho .* v_tip.^3 .* ss .* cd ./ dl ./ 8);
    else
        v_d = v_i * (k_i - 1.125 * v_y / v_i - 1.372 * (v_y / v_i)^2 - 1.718 * (v_y / v_i)^3 - 0.655 * (v_y / v_i)^4); % Induced velocity in descent according to an empirical relation (see lecture slides)
        pl = 1 ./ (v_y + sqrt(2) .* 2 .* k_i .* v_d + rho .* v_tip.^3 ./ dl .* ss .* cd ./ 8);
    end
else
    if v_y / v_i <= -2 % If this condition is met, the vertical climb equation is used for descent, else, an empirical equation is employed
        pl = 1 ./ (v_y - k_i ./ 2 * (v_y + sqrt(v_y.^2 - 2 .* dl ./ rho)) + rho .* v_tip.^3 .* ss .* cd ./ dl ./ 8);
    else
        v_d = v_i * (k_i - 1.125 * v_y / v_i - 1.372 * (v_y / v_i)^2 - 1.718 * (v_y / v_i)^3 - 0.655 * (v_y / v_i)^4); % Induced velocity in descent according to an empirical relation (see lecture slides)
        pl = 1 ./ (v_y + k_i .* v_d + rho .* v_tip.^3 ./ dl .* ss .* cd ./ 8);
    end
end

if strcmp(ducted, 'Yes')
    pl = pl .* sqrt(2);
end

% function pl = transition_constraint(wl, dl, rho, k, cd_0, v_tip, ss, cd, k_i, v, tt_tilt)
% aa = 0; % Assuming zero angle of attack of the blades
% mm = v * cosd(aa) / v_tip;
% pl = 1 ./ (k_i ./ sind(tt_tilt) .* sqrt(-v.^2 ./ 2 + sqrt((v.^2 ./ 2).^2 + (dl ./ 2 ./ rho ./ sind(tt_tilt)).^2)) + rho .* v_tip.^3 ./ dl .* (ss .* cd ./ 8 .* (1 + 4.6 .* mm.^2)) + 0.5 .* rho .* v^3 .* cd_0 ./ wl + 2 .* wl .* k ./ rho ./ v);

%% Vertical flight constraint regions
function c = hover_region(pl, dl, rho, m_i, mtom, v_tip, ss, cd, k_i, r, coaxial, ducted)
if strcmp(ducted, 'Yes')
    d = sqrt(2);
else
    d = 1;
end

if strcmp(coaxial, 'Yes')
    c = pl < d.* 1./ (1 ./ r .* m_i ./ mtom .* (k_i .* 2.* sqrt(2) .* sqrt(2 .* dl .* m_i ./ mtom ./ rho) + rho .* v_tip.^3 .* ss .* cd ./ (dl .* m_i ./ mtom .* 8)));
else
    c = pl < d.* 1./ (1 ./ r .* m_i ./ mtom .* (k_i .* sqrt(2 .* dl .* m_i ./ mtom ./ rho) + rho .* v_tip.^3 .* ss .* cd ./ (dl .* m_i ./ mtom .* 8)));
end

function c = vertical_climb_region(pl, dl, rho, m_i, mtom, v_tip, ss, cd, k_i, v_y, r, coaxial, ducted)
if strcmp(ducted, 'Yes')
    d = sqrt(2);
else
    d = 1;
end

if strcmp(coaxial, 'Yes')
    c = pl < d.* 1 ./ (1 ./ r .* m_i ./ mtom .*(v_y - 2.* sqrt(2) .* k_i .* v_y ./ 2 + 2.* sqrt(2) .* k_i .* sqrt(v_y.^2 + 2 .* dl .* m_i ./ mtom ./ rho) ./ 2 + rho .* v_tip.^3 .* ss .* cd ./ (dl .* m_i ./ mtom .* 8)));
else
    c = pl < d.* 1 ./ (1 ./ r .* m_i ./ mtom .*(v_y - k_i .* v_y ./ 2 + k_i .* sqrt(v_y.^2 + 2 .* dl .* m_i ./ mtom ./ rho) ./ 2 + rho .* v_tip.^3 .* ss .* cd ./ (dl .* m_i ./ mtom .* 8)));
end

function c = vertical_descent_region(pl, dl, rho, v_tip, ss, cd, k_i, v_y, v_i, coaxial, ducted)
if strcmp(coaxial, 'Yes')
    if v_y / v_i <= -2 % If this condition is met, the vertical climb equation is used for descent, else, an empirical equation is employed
        c = pl < 1 ./ (v_y - k_i ./ 2 * 2.* sqrt(2).* (v_y + sqrt(v_y.^2 - 2 .* dl ./ rho)) + rho .* v_tip.^3 .* ss .* cd ./ dl ./ 8);
    else
        v_d = v_i * (k_i - 1.125 * v_y / v_i - 1.372 * (v_y / v_i)^2 - 1.718 * (v_y / v_i)^3 - 0.655 * (v_y / v_i)^4); % Induced velocity in descent according to an empirical relation (see lecture slides)
        c = pl < 1 ./ (v_y + k_i .* 2.* sqrt(2).* v_d + rho .* v_tip.^3 ./ dl .* ss .* cd ./ 8);
    end
else
    if v_y / v_i <= -2 % If this condition is met, the vertical climb equation is used for descent, else, an empirical equation is employed
        c = pl < 1 ./ (v_y - k_i ./ 2 * 2.* sqrt(2) .* (v_y + sqrt(v_y.^2 - 2 .* dl ./ rho)) + rho .* v_tip.^3 .* ss .* cd ./ dl ./ 8);
    else
        v_d = v_i * (k_i - 1.125 * v_y / v_i - 1.372 * (v_y / v_i)^2 - 1.718 * (v_y / v_i)^3 - 0.655 * (v_y / v_i)^4); % Induced velocity in descent according to an empirical relation (see lecture slides)
        c = pl < 1 ./ (v_y + k_i .* 2.* sqrt(2) .* v_d + rho .* v_tip.^3 ./ dl .* ss .* cd ./ 8);
    end
end

if strcmp(ducted, 'Yes')
    c = sqrt(2) .* c;
end

% function c = transition_region(pl, wl, dl, rho, k, cd_0, v_tip, ss, cd, k_i, v, tt_tilt)
% aa = 0; % Assuming zero angle of attack of the blades
% mm = v * cosd(aa) / v_tip;
% c = pl < 1 ./ (k_i ./ sind(tt_tilt) .* sqrt(-v.^2 ./ 2 + sqrt((v.^2 ./ 2).^2 + (dl ./ 2 ./ rho ./ sind(tt_tilt)).^2)) + rho .* v_tip.^3 ./ dl .* (ss .* cd ./ 8 .* (1 + 4.6 .* mm.^2)) + 0.5 .* rho .* v^3 .* cd_0 ./ wl + 2 .* wl .* k ./ rho ./ v);

% Forward flight constraint functions
function wl = range_constraint_jet(rho, v, m_i, mtom, cd_0, k)
wl = mtom / m_i * 0.5 * rho * v^2 * sqrt(cd_0 / 3 / k);

function wl = range_constraint_prop(rho, v, m_i, mtom, cd_0, k)
wl = mtom / m_i * 0.5 * rho * v^2 * sqrt(cd_0 / k);

function wl = endurance_constraint_jet(rho, v, m_i, mtom, cd_0, k)
wl = mtom / m_i * 0.5 * rho * v^2 * sqrt(cd_0 / k);

function wl = endurance_constraint_prop(rho, v, m_i, mtom, cd_0, k)
wl =  mtom / m_i * 0.5 * rho * v^2 * sqrt(3 * cd_0 / k);

function wl = stall_speed_constraint(rho, v_s, m_i, mtom, c_lmax)
wl = mtom / m_i * 0.5 * rho * v_s^2 * c_lmax;

function ptl = cruise_speed_constraint_jet(wl, rho, v, m_i, mtom, cd_0, k)
ptl = 1 ./ (m_i ./ mtom .* (rho .* v.^2 .* cd_0 ./ 2 ./ wl .* m_i ./ mtom  + 2 .* k .* wl .* m_i ./ mtom ./ rho ./ v.^2));

function ptl = cruise_speed_constraint_prop(wl, rho, v, m_i, mtom, cd_0, k, ee, r)
ptl = ee .* r ./ (m_i ./ mtom .* (rho .* v.^3 .* cd_0 ./ 2 ./ wl .* m_i ./ mtom + 2 .* k .* wl .* m_i ./ mtom ./ rho ./ v));

function pl = fwd_flight_constraint_rotary_wing_prop(dl, a, s, v, r, ee, mr)
global constants
pl_i = 1 ./ (r{1}.induced_power_factor ./ (2 .* s.density .* s.velocity) .* (s.weight_init ./ v.mass).^2 .* dl);
pl_0 = 1./ (r{1}.rotor_solidity .* r{1}.base_drag_coefficient ./ 8 .* (1 + 4.65 .* advance_ratio(r{1}.rotational_speed,r{1}.radius,s.velocity,blade_aoa(0,a.lift_to_drag_ratio)).^2).*s.density ./ dl .* (r{1}.rotational_speed .* r{1}.radius).^3);
pl_p = 1./ (s.density(1) .* s.velocity.^3 .* a.equivalent_wetted_area ./ 2 ./ v.mass ./ constants.g);
pl = ee .* mr{1}.rating .* (pl_i + pl_0 + pl_p);

function pl = climb_constraint_jet(wl, rho, v, m_i, mtom, cd_0, k, gg)
pl = 1 ./ (m_i ./ mtom .*(sind(gg) + rho .* v.^2 .* cd_0 ./ 2 ./ wl .* m_i ./ mtom + 2 .* k .* wl * m_i ./ mtom ./ rho ./ v.^2));

function pl = climb_constraint_prop(wl, rho, v, m_i, mtom, cd_0, k, gg, ee, r)
pl = ee .* r ./ (m_i ./ mtom .* (v .* sind(gg) + rho .* v.^3 .* cd_0 ./ 2 ./ wl .* m_i ./ mtom + 2 .* k .* wl .* m_i ./ mtom ./ rho ./ v));

function pl = climb_constraint_rotary_prop(dl, a, s, v, r, ee, mr)
global constants
pl_i = 1 ./ (r{1}.induced_power_factor ./ (2 .* s.density(1) .* s.velocity) .* (s.weight_init ./ v.mass).^2 .* dl);
pl_0 = 1./ (r{1}.rotor_solidity .* r{1}.base_drag_coefficient ./ 8 .* (1 + 4.65 .* advance_ratio(r{1}.rotational_speed,r{1}.radius,s.velocity,blade_aoa(s.angle,a.lift_to_drag_ratio)).^2).*s.density(1) ./ dl .* (r{1}.rotational_speed .* r{1}.radius).^3);
pl_c = 1./ (s.velocity .* sind(s.angle) .* (s.weight_init ./ v.mass));
pl_p = 1./ (s.density(1) .* s.velocity.^3 .* a.equivalent_wetted_area ./ 2 ./ v.mass ./ constants.g);
pl = ee * mr{1}.rating * (pl_i + pl_0 + pl_c + pl_p);

function pl = climb_angle_constraint_jet(wl, rho, cd_0, k, gg)
pl = climb_constraint_jet(wl, rho, v_best_climb_angle_jet(wl, rho, k, cd_0), cd_0, k, gg); % TODO: Replace with segment speed

function pl = climb_angle_constraint_prop(wl, rho, cd_0, k, gg, ee)
pl = climb_constraint_prop(wl, rho, v_best_climb_angle_prop(wl, rho, k, cd_0), cd_0, k, gg, ee); % TODO: Replace with segment speed

% function pl = climb_rate(wl, rho, cd_0, k, gg, propulsion)
% if is_jet(propulsion)
%     % pl = fsolve(@(x)climb_rate_jet_error(x, wl, rho, cd_0, k, gg, propulsion), 0.01, optimoptions('fsolve', 'Display','iter'));
% elseif is_prop(propulsion)
%     pl = climb(wl, rho, v_best_climb_rate_prop(wl, rho, k, cd_0), cd_0, k, gg, propulsion);
% end

% function err = climb_rate_jet_error(tl, wl, rho, cd_0, k, gg, propulsion)
% err = climb(wl, rho, v_best_climb_rate_jet(tl, wl, rho, k, cd_0), cd_0, k, gg, propulsion) - tl;

% Forward flight constraint regions
function c = range_region_jet(wl, rho, v, m_i, mtom, cd_0, k)
c = wl < mtom / m_i * 0.5 * rho * v^2 * sqrt(cd_0 / 3 / k); 

function c = range_region_prop(wl, rho, v, m_i, mtom, cd_0, k)
c = wl < mtom / m_i * 0.5 * rho * v^2 * sqrt(cd_0 / k);

function c = endurance_region_jet(wl, rho, v, m_i, mtom, cd_0, k)
c = wl < mtom / m_i * 0.5 * rho * v^2 * sqrt(cd_0 / k);

function c = endurance_region_prop(wl, rho, v, m_i, mtom, cd_0, k)
c = wl < mtom / m_i * 0.5 * rho * v^2 * sqrt(3 * cd_0 / k);

function c = stall_speed_region(wl, rho, v_s, m_i, mtom, c_lmax)
c = wl < mtom / m_i * 0.5 * rho * v_s^2 * c_lmax;

function c = cruise_speed_region_jet(pl, wl, rho, v, m_i, mtom, cd_0, k)
c = pl < 1 ./ (m_i ./ mtom .* (rho .* v.^2 .* cd_0 ./ 2 ./ wl .* m_i ./ mtom  + 2 .* k .* wl .* m_i ./ mtom ./ rho ./ v.^2));

function c = cruise_speed_region_prop(pl, wl, rho, v, m_i, mtom, cd_0, k, ee, r)
c = pl < ee .* r ./ (m_i ./ mtom .* (rho .* v.^3 .* cd_0 ./ 2 ./ wl .* m_i ./ mtom + 2 .* k .* wl .* m_i ./ mtom ./ rho ./ v));

function c = fwd_flight_region_rotary_wing_prop(pl, dl, a, s, v, r, ee, mr)
global constants
pl_i = 1 ./ (r{1}.induced_power_factor ./ (2 .* s.density(1) .* s.velocity) .* (s.weight_init ./ v.mass).^2 .* dl);
pl_0 = 1./ (r{1}.rotor_solidity .* r{1}.base_drag_coefficient ./ 8 .* (1 + 4.65 .* advance_ratio(r{1}.rotational_speed,r{1}.radius,s.velocity,blade_aoa(0,a.lift_to_drag_ratio)).^2).*s.density(1) ./ dl .* (r{1}.rotational_speed .* r{1}.radius).^3);
pl_p = 1./ (s.density(1) .* s.velocity.^3 .* a.equivalent_wetted_area ./ 2 ./ v.mass ./ constants.g);
c = pl < ee .* mr{1}.rating .* (pl_i + pl_0 + pl_p);

function c = climb_region_jet(pl, wl, rho, v, m_i, mtom, cd_0, k, gg)
c = pl < 1 ./ (sind(gg) + rho .* v.^2 .* cd_0 ./ 2 ./ wl .* m_i ./ mtom + 2 .* k .* wl .* m_i ./ mtom ./ rho ./ v.^2);

function c = climb_region_prop(pl, wl, rho, v, m_i, mtom, cd_0, k, gg, ee, r)
c = pl < ee .* r ./ (v .* sind(gg) + rho .* v.^3 .* cd_0 ./ 2 ./ wl .* m_i ./ mtom + 2 .* k .* wl .* m_i ./ mtom ./ rho ./ v);

function c = climb_region_rotary_prop(pl, dl, a, s, v, r, ee, mr)
global constants
pl_i = 1 ./ (r{1}.induced_power_factor ./ (2 .* s.density(1) .* s.velocity) .* (s.weight_init ./ v.mass).^2 .* dl);
pl_0 = 1./ (r{1}.rotor_solidity .* r{1}.base_drag_coefficient ./ 8 .* (1 + 4.65 .* advance_ratio(r{1}.rotational_speed,r{1}.radius,s.velocity,blade_aoa(s.angle,a.lift_to_drag_ratio)).^2).*s.density(1) ./ dl .* (r{1}.rotational_speed .* r{1}.radius).^3);
pl_c = 1./ (s.velocity .* s.angle .* (s.weight_init ./ v.mass));
pl_p = 1./ (s.density(1) .* s.velocity.^3 .* a.equivalent_wetted_area ./ 2 ./ v.mass ./ constants.g);
c = pl < ee .* mr{1}.rating .* (pl_i + pl_0 + pl_c + pl_p);

function c = climb_angle_region_jet(pl, wl, rho, cd_0, k, gg)
c = climb_region_jet(pl, wl, rho, v_best_climb_angle_jet(wl, rho, k, cd_0), cd_0, k, gg);

function c = climb_angle_region_prop(pl, wl, rho, cd_0, k, gg, ee)
c = climb_region_prop(pl, wl, rho, v_best_climb_angle_prop(wl, rho, k, cd_0), cd_0, k, gg, ee);
