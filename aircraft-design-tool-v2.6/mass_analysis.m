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
       
function [mission, vehicle] = mass_analysis(mission, vehicle, energy)    

% Set energy source mass values to zero
for i = 1 : length(vehicle.components)
    if is_type(vehicle.components{i}, 'energy')
        vehicle.components{i}.mass = 0;
    end
end

mass_to = fsolve(@(x)mtow_error(x, mission, vehicle, energy), sum_masses(vehicle)); % iterative process to find mtom
[~, mission, vehicle] = mtow_error(mass_to, mission, vehicle, energy);
statistics(vehicle);

% Plot pie chart with mass distributions
function [] = statistics(vehicle)
payload = 0;
empty_mass = 0;
battery_mass = 0;
fuel_mass = 0;
hydrogen_mass = 0;

for i = 1 : length(vehicle.components)
    if is_type(vehicle.components{i},'mass.point')
        payload = payload + vehicle.components{i}.mass * vehicle.components{i}.number;
    elseif is_type(vehicle.components{i},'mass.empty') || is_type(vehicle.components{i},'fuselage') || is_type(vehicle.components{i},'wing') || is_type(vehicle.components{i},'landing_gear') 
        empty_mass = empty_mass + vehicle.components{i}.mass;
    elseif is_type(vehicle.components{i},'boom') || is_type(vehicle.components{i},'tank') || is_type(vehicle.components{i},'engine') || is_type(vehicle.components{i},'driver') || is_type(vehicle.components{i},'generator') || is_type(vehicle.components{i},'gearbox') || is_type(vehicle.components{i},'motor') || is_type(vehicle.components{i},'fuel_cell') || is_type(vehicle.components{i},'valve') || is_type(vehicle.components{i},'inverter') || is_type(vehicle.components{i},'power_dist_module') || is_type(vehicle.components{i},'electronic_speed_controller')
        empty_mass = empty_mass + vehicle.components{i}.mass * vehicle.components{i}.number;
    elseif is_type(vehicle.components{i},'energy.fuel')
        fuel_mass = fuel_mass + vehicle.components{i}.mass;
    elseif is_type(vehicle.components{i},'energy.electric')
        battery_mass = battery_mass + vehicle.components{i}.mass;
    elseif is_type(vehicle.components{i},'energy.hydrogen')
        hydrogen_mass = hydrogen_mass + vehicle.components{i}.mass;
    end
end

labels = {'Payload','Empty','Battery','Fuel','Hydrogen'};
masses = [payload,empty_mass,battery_mass,fuel_mass,hydrogen_mass];

figure
pie(masses)
legend(labels,'Location','bestoutside','Orientation','horizontal')

% Sum masses of all components
function mass = sum_masses(vehicle)
mass = 0;
for i = 1 : length(vehicle.components)
    m = vehicle.components{i}.mass;

    if isfield(vehicle.components{i}, 'number')
        m = m * vehicle.components{i}.number;
    end
       
    mass = mass + m;
end

% Breguet equation for all cruise segments
function mf_fuel = breguet(range, velocity, sfc, ld, fuel_weight_fraction, battery_energy_fraction)
global constants;
mf_fuel = 1 - exp(-range * sfc * constants.g * (1 - battery_energy_fraction)/ (velocity * ld * (1 - fuel_weight_fraction)));

function dyn_press = dynamic_pressure(segment)
dyn_press = 0.5 * segment.density * segment.velocity^2;

% Lift-to-drag ratio of aircraft
function ld = get_ld(vehicle, segment, segment_id)
global constants;

w = find_by_type(vehicle.components, 'wing.main');
if isfield(w{1},'lift_to_drag_ratio')
    ld = w{1}.lift_to_drag_ratio;
else
parasitic_drag_wing = dynamic_pressure(segment) * w{1}.area_ref * w{1}.segments{segment_id}.base_drag_coefficient;
for i = 1 : length(vehicle.components)
    if is_type(vehicle.components{i}, 'wing.htail') || is_type(vehicle.components{i}, 'wing.canard') || is_type(vehicle.components{i}, 'wing.secondary')
        s = vehicle.components{i};
        parasitic_drag_stabilizer = dynamic_pressure(segment) * s.area_ref * s.segments{segment_id}.base_drag_coefficient;
    end
end
mu = s.span/w{1}.span;
total_drag = parasitic_drag_wing + parasitic_drag_stabilizer + (segment.weight_init * constants.g)^2/(pi() * dynamic_pressure(segment) * w{1}.span^2)*((1-w{1}.mutual_interf^2)/(1-2*w{1}.mutual_interf*mu+mu^2));
ld = segment.weight_init * constants.g / total_drag;
end

% Taxi segment calculations - Fixed-wing
function [vehicle, segment] = taxi(segment, vehicle, energy)
global databases

[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');

segment.weight_init = vehicle.mass;
if is_type(source, 'energy.fuel')
    jet_fuel_id = find_data('jet-fuel','source');
    fuel_type_id = find_data(source{1}.source_type,'source');
    mf_fuel = (1 - 0.9725)*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    vehicle.components{network_ids(source_id)}.mass = source{1}.mass + mf_fuel * vehicle.mass;
    vehicle.mass = vehicle.mass * (1 - mf_fuel);
end
segment.weight_final = vehicle.mass;

% Hover segment calculations - VTOL
function [vehicle, segment, energy_required,energy_current] = hover(segment, vehicle, energy, energy_required,energy_current)
global databases
segment.weight_init = vehicle.mass;

[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
rotor = find_by_type(network, 'driver.prop');

% Network efficiency
e = network_efficiency(network, 'no_type');

% Calculate thrust required per rotor
rotor = balance_forces_vertical(rotor, vehicle.mass);

battery_mass = 0;

[source, source_id] = find_by_type(network, 'energy');
if is_type(source, 'energy.fuel')
    errordlg('Hover not available for fuel energy sources'); % NOT AVAILABLE
    return;
elseif is_type(source, 'energy.electric')
    battery_type_id = find_data(source{1}.source_type,'source');
    for i = 1 : size(rotor,1)
        induced_power = rotor{i,1}.induced_power_factor * rotor{i,1}.number * rotor{i,1}.thrust * sqrt(rotor{i,1}.thrust / (2 * segment.density * rotor_area(rotor{i,1}))); % total induced power
        profile_power = rotor_area(rotor{i,1}) * segment.density * rotor{i,1}.number * rotor{i,1}.tip_velocity^3 * rotor{i,1}.rotor_solidity * rotor{i,1}.base_drag_coefficient / 8; % profile power for each rotor
        if strcmp(rotor{i,1}.coaxial,'Yes')
            induced_power = induced_power * 2*sqrt(2);
            profile_power = profile_power * 2;
        end
        if strcmp(rotor{i,1}.ducted,'Yes')
            induced_power = induced_power/sqrt(2);
        end
        battery_mass = battery_mass + (induced_power + profile_power) * segment.time / databases.source.specific_energy{battery_type_id,1} / e(i,1) / source{1}.usable_fraction;
    end
    energy_current = energy_current + battery_mass * databases.source.specific_energy{battery_type_id,1}; 
    if energy_current >= energy_required
        vehicle.components{network_ids(source_id)}.mass = source{1}.mass + (energy_current - energy_required) / databases.source.specific_energy{battery_type_id,1};
        energy_required = energy_current;
    elseif energy_current < energy_required
        % Do not add any battery mass
    end
end
   
segment.weight_final = vehicle.mass;

% Climb segment calculations - Fixed-wing and rotary-wing
function [vehicle, segment, energy_required, energy_current] = climb(segment, vehicle, energy, energy_required, energy_current)
global constants databases
[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');
[segment_props, ~] = find_by_name(vehicle.segments, segment.name);
aircraft = find_by_type(vehicle.components, 'aircraft');
rotor = find_by_type(network, 'driver.prop');

segment.weight_init = vehicle.mass;
if is_type(source, 'energy.fuel')
    mach = segment.velocity / segment.speed_sound(2);
    jet_fuel_id = find_data('jet-fuel','source');
    fuel_type_id = find_data(source{1}.source_type,'source');
    if mach < 1
        mf_fuel = (1 - (1 - 0.04 * mach))*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    else
        mf_fuel = (1 - (0.96 - 0.03 * (mach - 1)))*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    end
    vehicle.components{network_ids(source_id)}.mass = source{1}.mass + mf_fuel * vehicle.mass;
    vehicle.mass = vehicle.mass - mf_fuel * vehicle.mass;
elseif is_type(source, 'energy.electric')
    battery_type_id = find_data(source{1}.source_type,'source');
    if is_type(aircraft, 'aircraft.rotary_wing') % Climb segment for rotary-wing aircraft (TEMP: assuming all rotors are of the same size)
        induced_power = rotor{1}.induced_power_factor * (vehicle.mass * constants.g)^2 / 2 / segment.density(1) / rotor{1}.number / rotor_area(rotor{1}) / segment.velocity;
        profile_power = rotor_area(rotor{1}) * rotor{1}.number * segment.density(1) * rotor{1}.rotor_solidity * rotor{1}.base_drag_coefficient * (1 + 4.65 * advance_ratio(rotor{1}.rotational_speed, rotor{1}.radius, segment.velocity, blade_aoa(segment.angle, aircraft{1}.rotary_wing.lift_to_drag_ratio))^2) * (rotor{1}.rotational_speed * rotor{1}.radius)^3 / 8;
        climb_power = vehicle.mass * constants.g * segment.velocity * sind(segment.angle);
        parasitic_power = 0.5 * segment.density(1) * segment.velocity^3 * aircraft{1}.rotary_wing.equivalent_wetted_area;
        total_power_climb = induced_power + profile_power + climb_power + parasitic_power;
        m_batt = total_power_climb * segment.time / network_efficiency(network, 'no_type') / source{1}.usable_fraction / databases.source.specific_energy{battery_type_id,1};
        energy_current = energy_current + m_batt * databases.source.specific_energy{battery_type_id,1};

    else % Climb segment for fixed-wing aircraft
        term_2 = segment.density(1) * segment.velocity^2 * segment_props.base_drag_coefficient / (2 * (vehicle.mass * constants.g / vehicle.area_ref));
        c = find_by_type(vehicle.components, 'wing.main');
        k = 1 / pi / c{1}.aspect_ratio / c{1}.oswald_efficiency;
        term_3 = 2 * k / segment.density(1) / segment.velocity^2 * vehicle.mass * constants.g / vehicle.area_ref;
        p_w = segment.velocity * (sind(segment.angle) + term_2 + term_3);
        mf_batt = p_w * constants.g * (segment.altitude(2) - segment.altitude(1)) / databases.source.specific_energy{battery_type_id,1} / network_efficiency(network, 'no_type') / segment.velocity / sind(segment.angle) / source{1}.usable_fraction;
        energy_current = energy_current + mf_batt * vehicle.mass * databases.source.specific_energy{battery_type_id,1};
    end

    if energy_current >= energy_required
        vehicle.components{network_ids(source_id)}.mass = source{1,1}.mass + (energy_current - energy_required) / databases.source.specific_energy{battery_type_id,1};
        energy_required = energy_current;
    elseif energy_current < energy_required
        % Do not add any battery mass
    end
end
segment.weight_final = vehicle.mass;

% Vertical climb segment calculations - VTOL
function [vehicle, segment, energy_required, energy_current] = vertical_climb(segment, vehicle, energy, energy_required, energy_current)
global databases 
segment.weight_init = vehicle.mass;

[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
rotor = find_by_type(network, 'driver.prop');

% Network efficiency
e = network_efficiency(network, 'no_type'); % must include all components, including batteries

% Calculate thrust required per rotor
rotor = balance_forces_vertical(rotor, vehicle.mass);

battery_mass = 0;

[source, source_id] = find_by_type(network, 'energy');
if is_type(source, 'energy.fuel')
    errordlg('Vertical climb not available for fuel energy sources'); % NOT AVAILABLE
    return;
elseif is_type(source, 'energy.electric')
    battery_type_id = find_data(source{1}.source_type,'source');
    altitude_range = segment.altitude(2) - segment.altitude(1);
    for i = 1 : size(rotor,1)
        induced_velocity = -0.5 * segment.velocity + sqrt((segment.velocity/2)^2 + rotor{i,1}.thrust / (2 * segment.density(1) * rotor_area(rotor{i,1})));
        if strcmp(rotor{i,1}.coaxial,'Yes')
            induced_velocity = induced_velocity * 2*sqrt(2);
        end
        if strcmp(rotor{i,1}.ducted,'Yes')
            induced_velocity = induced_velocity/sqrt(2);
        end
        induced_power = rotor{i,1}.number * rotor{i,1}.thrust * (segment.velocity + rotor{i,1}.induced_power_factor * induced_velocity);
        profile_power = rotor_area(rotor{i,1}) * segment.density(1) * rotor{i,1}.number * rotor{i,1}.tip_velocity^3 * rotor{i,1}.rotor_solidity * rotor{i,1}.base_drag_coefficient / 8; 
        if strcmp(rotor{i,1}.coaxial,'Yes')
            profile_power = profile_power * 2;
        end
        battery_mass = battery_mass + (induced_power + profile_power) * altitude_range / databases.source.specific_energy{battery_type_id,1} / e(i,1) / source{1,1}.usable_fraction / segment.velocity;
    end

    energy_current = energy_current + battery_mass * databases.source.specific_energy{battery_type_id,1};
    if energy_current >= energy_required
        vehicle.components{network_ids(source_id)}.mass = source{1,1}.mass + (energy_current - energy_required) / databases.source.specific_energy{battery_type_id,1};
        energy_required = energy_current;
    elseif energy_current < energy_required
        % Do not add any battery mass
    end
end
  
segment.weight_final = vehicle.mass;

% Acceleration segment calculations - Fixed-wing
function [vehicle, segment, energy_required,energy_current] = acceleration(segment, prev_segment, vehicle, energy, energy_required,energy_current)
global databases
[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');
segment.weight_init = vehicle.mass;
mach = segment.velocity / segment.speed_sound;

if is_type(source, 'energy.fuel')
    jet_fuel_id = find_data('jet-fuel','source');
    fuel_type_id = find_data(source{1}.source_type,'source');
    if mach == prev_segment.velocity / prev_segment.speed_sound
        mf_fuel = 0;
    elseif mach < 1
        mf_fuel = (1 - (1 - 0.04 * mach))*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    else
        mf_fuel = (1 - (0.96 - 0.03 * (mach - 1)))*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    end
    vehicle.components{network_ids(source_id)}.mass = source.mass + mf_fuel * vehicle.mass;
    vehicle.mass = vehicle.mass - mf_fuel * vehicle.mass;
elseif is_type(source, 'energy.electric')
    errordlg('Acceleration not available for electric energy sources'); % NOT AVAILABLE
    return;
end
segment.weight_final = vehicle.mass;

% Cruise segment calculations - Fixed-wing
function [vehicle, segment, energy_required,energy_current] = cruise(segment, segment_id, vehicle, energy, energy_required,energy_current)
global constants databases;
segment.weight_init = vehicle.mass;
[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');
aircraft = find_by_type(vehicle.components, 'aircraft');
rotor = find_by_type(network, 'driver.prop');

if is_type(source, 'energy.fuel')
    ld = get_ld(vehicle, segment, segment_id);
    jet_fuel_id = find_data('jet-fuel','source');
    fuel_type_id = find_data(source{1}.source_type,'source');
    engine = find_by_type(network, 'engine');
    if is_type(engine, 'engine.jet')
        mf_fuel = breguet(segment.range, segment.velocity, engine{1}.specific_fuel_consumption, ld,0,0)*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    elseif is_type(engine, 'engine.prop')
        equivalent_sfc = engine{1,1}.brake_specific_fuel_consumption * segment.velocity / network_efficiency(network, 'no_type');
        mf_fuel = breguet(segment.range, segment.velocity, equivalent_sfc, ld,0,0) * databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    end
    vehicle.components{network_ids(source_id)}.mass = source{1,1}.mass + mf_fuel * vehicle.mass;
    vehicle.mass = vehicle.mass - mf_fuel * vehicle.mass;
elseif is_type(source, 'energy.electric')
    battery_type_id = find_data(source{1}.source_type,'source');
    if is_type(aircraft, 'aircraft.rotary_wing') % Cruise segment for rotary-wing aircraft (TEMP: assuming all rotors are of the same size)
        induced_power = rotor{1}.induced_power_factor * (vehicle.mass * constants.g)^2 / 2 / segment.density / rotor{1}.number / rotor_area(rotor{1}) / segment.velocity;
        profile_power = rotor_area(rotor{1}) * rotor{1}.number * segment.density * rotor{1}.rotor_solidity * rotor{1}.base_drag_coefficient * (1 + 4.65 * advance_ratio(rotor{1}.rotational_speed, rotor{1}.radius, segment.velocity, blade_aoa(0, aircraft{1}.rotary_wing.lift_to_drag_ratio))^2) * (rotor{1}.rotational_speed * rotor{1}.radius)^3 / 8;
        parasitic_power = 0.5 * segment.density * segment.velocity^3 * aircraft{1}.rotary_wing.equivalent_wetted_area;
        total_power_cruise = induced_power + profile_power + parasitic_power;
        m_batt = total_power_cruise * segment.time / network_efficiency(network, 'no_type') / source{1}.usable_fraction / databases.source.specific_energy{battery_type_id,1};
        energy_current = energy_current + m_batt * databases.source.specific_energy{battery_type_id,1};
    else
        ld = get_ld(vehicle, segment, segment_id);
        mf_batt = segment.range * constants.g / databases.source.specific_energy{battery_type_id,1} / network_efficiency(network, 'no_type') / ld / source{1}.usable_fraction;
        energy_current = energy_current + mf_batt * vehicle.mass * databases.source.specific_energy{battery_type_id,1};
    end    

    if energy_current >= energy_required
        vehicle.components{network_ids(source_id)}.mass = source{1}.mass + (energy_current - energy_required) / databases.source.specific_energy{battery_type_id,1};
        energy_required = energy_current;
    elseif energy_current < energy_required
        % Do not add any battery mass
    end
end
segment.weight_final = vehicle.mass;

% Hybrid cruise segment calculations - Fixed-wing
function [vehicle, segment, energy_required, energy_current] = hybrid_cruise(segment, segment_id, vehicle, energy, energy_required,energy_current)
global constants databases;
segment.weight_init = vehicle.mass;

[battery_network, battery_network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network_battery));
[battery_source, battery_source_id] = find_by_type(battery_network, 'energy');

[fuel_network, fuel_network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network_fuel));
[fuel_source, fuel_source_id] = find_by_type(fuel_network, 'energy');

ld = get_ld(vehicle, segment, segment_id);
fuel_type_id = find_data(fuel_source{1}.source_type,'source');
battery_type_id = find_data(battery_source{1}.source_type,'source');

if is_type(fuel_source,'energy.fuel') % for jet-fuel/SAF 
    jet_fuel_id = find_data('jet-fuel','source');

    engine = find_by_type(fuel_network, 'engine');

    if is_type(engine, 'engine.jet')
        errordlg('Jet engine not available for recharging segments.'); % NOT AVAILABLE
    elseif is_type(engine, 'engine.prop')
        equivalent_sfc = engine{1}.brake_specific_fuel_consumption * segment.velocity / network_efficiency(fuel_network, 'no_type');
        mf_fuel = breguet(segment.range, segment.velocity, equivalent_sfc, ld, 0, battery_source{1}.battery_energy_fraction)*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    end
    vehicle.components{fuel_network_ids(fuel_source_id)}.mass = fuel_source{1}.mass + mf_fuel * vehicle.mass;
    term_1 = (battery_source{1}.battery_energy_fraction * segment.velocity)/ network_efficiency(battery_network, 'no_type') / battery_source{1}.usable_fraction;
    term_2 = 1/(equivalent_sfc * (1 - battery_source{1}.battery_energy_fraction));
    term_3 = exp(-segment.range * equivalent_sfc  * constants.g * (1 - battery_source{1}.battery_energy_fraction) / (segment.velocity * ld ));
    mf_batt = (term_1 * term_2 * (1- term_3)) / databases.source.specific_energy{battery_type_id,1};
    
    vehicle.mass = vehicle.mass - mf_fuel * vehicle.mass;
    segment.weight_final = vehicle.mass;
elseif is_type(fuel_source,'energy.hydrogen') % for fuel cell systems
    vehicle.components{fuel_network_ids(fuel_source_id)}.mass = vehicle.mass * constants.g * segment.range*(1-battery_source{1}.battery_energy_fraction)/(ld * network_efficiency(fuel_network, 'no_type') * databases.source.specific_energy{fuel_type_id,1}); % adjust variables
    mf_batt = constants.g * segment.range*battery_source{1}.battery_energy_fraction/(ld * network_efficiency(battery_network, 'no_type') * databases.source.specific_energy{battery_type_id,1}); % adjust variable
end

energy_current = energy_current + mf_batt * vehicle.mass * databases.source.specific_energy{battery_type_id,1};
if energy_current >= energy_required
    vehicle.components{battery_network_ids(battery_source_id)}.mass = battery_source{1}.mass + (energy_current - energy_required) / databases.source.specific_energy{battery_type_id,1};
    energy_required = energy_current;
elseif energy_current < energy_required
    %Do not add any battery mass
end

% Battery recharge cruise segment calculations - Fixed-wing
function [vehicle, segment, energy_required,energy_current] = recharge_cruise(segment, segment_id, vehicle, energy, energy_required,energy_current)
global databases;
segment.weight_init = vehicle.mass;
[prop_network, prop_network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network_prop));
[source, source_id] = find_by_type(prop_network, 'energy');

[recharge_network, ~] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network_recharge));

engine = find_by_type(prop_network, 'engine');
ld = get_ld(vehicle, segment, segment_id);
engine_recharge = find_by_type(recharge_network, 'engine');

if is_type(source, 'energy.fuel')
    jet_fuel_id = find_data('jet-fuel','source');
    fuel_type_id = find_data(source{1}.source_type,'source');
    if is_type(engine, 'engine.jet')
        errordlg('Jet engine not available for recharging segments.'); % NOT AVAILABLE
    elseif is_type(engine, 'engine.prop')
        equivalent_sfc = engine{1}.brake_specific_fuel_consumption * segment.velocity / network_efficiency(prop_network, 'no_type');
        mf_fuel = breguet(segment.range, segment.velocity, equivalent_sfc, ld, engine_recharge{1,1}.fuel_weight_fraction,0)*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    end
    vehicle.components{prop_network_ids(source_id)}.mass = source{1}.mass + mf_fuel * vehicle.mass;
    
    energy_fuel = engine_recharge{1}.fuel_weight_fraction* mf_fuel * vehicle.mass * databases.source.specific_energy{fuel_type_id,1} * network_efficiency(recharge_network, 'no_type');
    
    if energy_fuel <= energy_current
        energy_current = energy_current - energy_fuel;
    else
        energy_current = 0;
        % errordlg('ATTENTION: You are exceeding the amount of energy needed to fully recharge your batteries! Consider shortening the duration of this flight segment.');
    end
    
    vehicle.mass = vehicle.mass - mf_fuel * vehicle.mass;
    
elseif is_type(source, 'energy.electric')
     errordlg('You must define a fuel energy source to recharge the batteries.'); % NOT AVAILABLE
     return;
end
segment.weight_final = vehicle.mass;

% Hold/endurance/loiter segment calculations - Fixed-wing
function [vehicle, segment, energy_required,energy_current] = hold(segment, segment_id, vehicle, energy, energy_required,energy_current)
global constants databases;
segment.weight_init = vehicle.mass;
[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');
engine = find_by_type(network, 'engine');

ld = get_ld(vehicle, segment, segment_id);

if is_type(source, 'energy.fuel')
    jet_fuel_id = find_data('jet-fuel','source');
    fuel_type_id = find_data(source{1}.source_type,'source');
    if is_type(engine, 'engine.jet')
        mf_fuel = breguet(segment.range, segment.velocity, engine{1}.specific_fuel_consumption, ld,0,0)*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    elseif is_type(engine, 'engine.prop')
        equivalent_sfc = engine{1}.brake_specific_fuel_consumption * segment.velocity / network_efficiency(network,'no_type');
        mf_fuel = breguet(segment.range, segment.velocity, equivalent_sfc, ld,0,0)*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    end
    vehicle.components{network_ids(source_id)}.mass = source{1}.mass + mf_fuel * vehicle.mass;
    vehicle.mass = vehicle.mass - mf_fuel * vehicle.mass;
elseif is_type(source, 'energy.electric')
    battery_type_id = find_data(source{1}.source_type,'source');
    mf_batt = segment.time * segment.velocity * constants.g / databases.source.specific_energy{battery_type_id,1} / network_efficiency(network,'no_type') / ld / source{1}.usable_fraction;
        
    energy_current = energy_current + mf_batt * vehicle.mass * databases.source.specific_energy{battery_type_id,1};
    if energy_current >= energy_required
        vehicle.components{network_ids(source_id)}.mass = source{1}.mass + (energy_current - energy_required) / databases.source.specific_energy{battery_type_id,1};
        energy_required = energy_current;
    elseif energy_current < energy_required
        % Do not add any battery mass
    end
end
segment.weight_final = vehicle.mass;

% Descent (glide) segment calculations - Fixed-wing and rotary-wing
% aircraft
% TO DO: determine descent angle for rotary-wing aircraft using iterative process
function [vehicle, segment] = descent(segment, vehicle)
segment.weight_init = vehicle.mass;
segment.weight_final = vehicle.mass;

% Vertical descent segment calculations - VTOL
function [vehicle, segment, energy_required,energy_current] = vertical_descent(segment, vehicle, energy, energy_required,energy_current)
global databases 
segment.weight_initial = vehicle.mass;

[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
rotor = find_by_type(network, 'driver.prop');

% Calculate thrust required per rotor
rotor = balance_forces_vertical(rotor, vehicle.mass);

% Network efficiency
e = network_efficiency(network, 'no_type');

battery_mass = 0;

[source, source_id] = find_by_type(network, 'energy');
if is_type(source, 'energy.fuel')
    errordlg('Vertical descent not available for fuel energy sources'); % NOT AVAILABLE
    return;
elseif is_type(source, 'energy.electric')
    battery_type_id = find_data(source{1}.source_type,'source');
    altitude_range = abs(segment.altitude(2) - segment.altitude(1));
    for i = 1 : size(rotor,1)
        v_i = sqrt(rotor{i,1}.thrust / 2 / segment.density(2) / rotor_area(rotor{i,1}) / rotor{i,1}.number); % Induced velocity in hover
        if segment.velocity / v_i <= -2 % If this condition is met, the vertical climb equation is used for descent, else, an empirical equation is employed
            v_d = -0.5 * segment.velocity - sqrt((0.5 * segment.velocity)^2 - rotor{i,1}.thrust / segment.density(2) / 2);
        else 
            v_d = v_i * (rotor{i,1}.induced_power_factor - 1.125 * segment.velocity / v_i - 1.372 * (segment.velocity / v_i)^2 - 1.718 * (segment.velocity / v_i)^3 - 0.655 * (segment.velocity / v_i)^4); % Induced velocity in descent according to an empirical relation (see lecture slides)
        end
        
        if strcmp(rotor{i,1}.coaxial,'Yes') % if rotors are coaxial, affect induced power through induced velocity
            v_d = v_d * 2*sqrt(2);
        end
        if strcmp(rotor{i,1}.ducted,'Yes') % if rotors are ducted, affect induced power through induced velocity
            v_d = v_d/sqrt(2);
        end
        induced_power = rotor{i,1}.number * rotor{i,1}.thrust * (segment.velocity + rotor{i,1}.induced_power_factor * v_d);
        profile_power = rotor_area(rotor{i,1}) * segment.density(2) * rotor{i,1}.number * rotor{i,1}.tip_velocity^3 * rotor{i,1}.rotor_solidity * rotor{i,1}.base_drag_coefficient / 8; 

        if strcmp(rotor{i,1}.coaxial,'Yes') % if rotors are coaxial, profile power is doubled
            profile_power = profile_power * 2;
        end

        if (induced_power + profile_power) > 0
            battery_mass = battery_mass + (induced_power + profile_power) * altitude_range / databases.source.specific_energy{battery_type_id,1} / e(i,1) / source{1}.usable_fraction / abs(segment.velocity);
        else
            battery_mass = 0;
        end
    end
    
    energy_current = energy_current + battery_mass * databases.source.specific_energy{battery_type_id,1};
    if energy_current >= energy_required
        vehicle.components{network_ids(source_id)}.mass = source{1}.mass + (energy_current - energy_required) / databases.source.specific_energy{battery_type_id,1};
        energy_required = energy_current;
    elseif energy_current < energy_required
        % Do not add any battery mass
    end
end  

segment.weight_final = vehicle.mass;

% Landing segment calculations - Fixed-wing
function [vehicle, segment] = landing(segment, vehicle, energy)
global databases
[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');
segment.weight_init = vehicle.mass;
if is_type(source, 'energy.fuel')
    jet_fuel_id = find_data('jet-fuel','source');
    fuel_type_id = find_data(source{1}.source_type,'source');
    mf_fuel = (1 - 0.9725)*databases.source.specific_energy{jet_fuel_id,1}/databases.source.specific_energy{fuel_type_id,1};
    vehicle.components{network_ids(source_id)}.mass = source{1}.mass + mf_fuel * vehicle.mass;
    vehicle.mass = vehicle.mass - mf_fuel * vehicle.mass;
elseif is_type(source, 'energy.electric')
    % mf_batt = 0;
end
segment.weight_final = vehicle.mass;

% Load/drop segment calculations - Fixed-wing
function [vehicle, segment] = load_step(vehicle, segment)
segment.weight_init = vehicle.mass;
vehicle.mass = vehicle.mass + segment.mass;
segment.weight_final = vehicle.mass;

% Function that computes difference between initial and final mtow
% estimates
function [error, mission, vehicle] = mtow_error(x, mission, vehicle, energy)
global databases 
vehicle.mass = x;

% Calculate empty weight
[empty, empty_id] = find_by_type(vehicle.components, 'mass.empty');
if (~isempty(empty))
    k = 1;
    if strcmp(empty{1}.material, 'composites') % if structure is mostly using composite materials
        k = 0.95;
    end
    if strcmp(empty{1}.equation.type, 'ax^(1+c)') % power law
        vehicle.components{empty_id}.mass = k * empty{1}.equation.a * vehicle.mass^(empty{1}.equation.c + 1);
    elseif strcmp(empty{1}.equation.type, 'ax+c') % linear
        vehicle.components{empty_id}.mass = k * empty{1}.equation.a * vehicle.mass + empty{1}.equation.c;
    end
end

% Initialize required and current battery energies to zero
energy_required = 0;
energy_current = 0;

% Iterate over mission segments
for i = 1 : length(mission.segments)
    [mission.segments{i}.temperature, mission.segments{i}.speed_sound, mission.segments{i}.pressure, mission.segments{i}.density] = atmosisa(mission.segments{i}.altitude);

    if strcmp(mission.segments{i}.type, 'taxi') % Taxi segment
        [vehicle, mission.segments{i}] = taxi(mission.segments{i}, vehicle, energy);
    elseif strcmp(mission.segments{i}.type, 'hover') % Hover segment
        [vehicle, mission.segments{i}, energy_required,energy_current] = hover(mission.segments{i}, vehicle, energy, energy_required,energy_current);
    elseif strcmp(mission.segments{i}.type, 'climb') % Climb segment
        [vehicle, mission.segments{i}, energy_required,energy_current]= climb(mission.segments{i}, vehicle, energy, energy_required,energy_current);
    elseif strcmp(mission.segments{i}.type, 'vertical_climb') % Vertical climb segment
        [vehicle, mission.segments{i}, energy_required,energy_current] = vertical_climb(mission.segments{i}, vehicle, energy, energy_required,energy_current);
    elseif strcmp(mission.segments{i}.type, 'acceleration') % Acceleration segment
        [vehicle, mission.segments{i} ,energy_required,energy_current] = acceleration(mission.segments{i}, mission.segments{i - 1}, vehicle, energy, energy_required,energy_current);
    elseif strcmp(mission.segments{i}.type, 'cruise') % Cruise segment
        [vehicle, mission.segments{i},energy_required,energy_current] = cruise(mission.segments{i}, i, vehicle, energy, energy_required,energy_current);
    elseif strcmp(mission.segments{i}.type, 'recharge_cruise') % Battery recharging cruise segment
        [vehicle, mission.segments{i},energy_required,energy_current] = recharge_cruise(mission.segments{i}, i, vehicle, energy, energy_required,energy_current);
    elseif strcmp(mission.segments{i}.type, 'hybrid_cruise') % Hybrid segment using fuel and battery as energy sources
        [vehicle, mission.segments{i},energy_required,energy_current] = hybrid_cruise(mission.segments{i}, i, vehicle, energy, energy_required,energy_current);
    elseif strcmp(mission.segments{i}.type, 'hold') % Hold segment
        [vehicle, mission.segments{i},energy_required,energy_current] = hold(mission.segments{i}, i, vehicle, energy, energy_required,energy_current);
    elseif strcmp(mission.segments{i}.type, 'combat') % Combat segment
    %     ei = find_energy_source_index(vehicle, mission.segments{i}.source);
    %     if is_type(source, 'energy.fuel')
    %         % mf_fuel = 1 - % TODO
    %     elseif is_type(source, 'energy.electric')
    %         errordlg('Combat not available for electric energy sources'); % NOT AVAILABLE
    %         break;
    %     end
    elseif strcmp(mission.segments{i}.type, 'descent') % Descent segment
        [vehicle, mission.segments{i}] = descent(mission.segments{i}, vehicle);
    elseif strcmp(mission.segments{i}.type, 'vertical_descent') % Vertical descent segment
        [vehicle, mission.segments{i},energy_required,energy_current] = vertical_descent(mission.segments{i}, vehicle, energy, energy_required,energy_current);
    elseif strcmp(mission.segments{i}.type, 'landing') % Landing segment
        [vehicle, mission.segments{i}] = landing(mission.segments{i}, vehicle, energy);
    elseif strcmp(mission.segments{i}.type, 'load_step') % Load step segment
        [vehicle, mission.segments{i}] = load_step(vehicle, mission.segments{i});
    end
end

for i = 1 : length(vehicle.components)
    if isfield(vehicle.components{i}, 'reserve')
        vehicle.components{i}.reserve_mass = vehicle.components{i}.mass * vehicle.components{i}.reserve;
        vehicle.components{i}.mass = vehicle.components{i}.mass + vehicle.components{i}.reserve_mass;
    end
    
    if is_type(vehicle.components{i}, 'energy.electric')
        battery_type_id = find_data(vehicle.components{i}.source_type,'source');
        reserve_energy = vehicle.components{i}.reserve_mass * databases.source.specific_energy{battery_type_id,1};
        total_energy = vehicle.components{i}.mass * databases.source.specific_energy{battery_type_id,1};
        vehicle.components{i}.remaining_energy = energy_required - energy_current + reserve_energy; 
        vehicle.components{i}.final_state_of_charge = vehicle.components{i}.remaining_energy / total_energy * 100; % we simplify by using the ratio of energies
        vehicle.components{i}.cell_capacity = (total_energy /(vehicle.components{i}.volt * vehicle.components{i}.number_packs.series * 3600))/vehicle.components{i}.number_packs.parallel; % Calculate amp-hour per battery cell
    end
end

vehicle.mass = x;

% Accumulate component masses and calculate error
error = vehicle.mass - sum_masses(vehicle);
