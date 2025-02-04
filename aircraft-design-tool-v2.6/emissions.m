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

function [vehicle,mission] = emissions(mission, vehicle)
global databases

[fuel, fuel_id] = find_by_type(vehicle.components, 'energy.fuel');
[hydrogen, hydrogen_id] = find_by_type(vehicle.components, 'energy.hydrogen');
[electric, electric_id] = find_by_type(vehicle.components, 'energy.electric');

if ~isempty(fuel) % emissions during operation and production
    jet_fuel_id = find_data(fuel{1,1}.source_type,'source');
    vehicle.components{fuel_id}.emissions.operation = (fuel{1}.mass-fuel{1}.reserve_mass)*databases.source.operation_EI_CO2{jet_fuel_id};
    vehicle.components{fuel_id}.emissions.production = fuel{1}.mass*databases.source.production_EI_CO2_per_kwh{jet_fuel_id}*databases.source.specific_energy{jet_fuel_id}/(3.6e6);
end

if ~isempty(electric) % emissions during operation and production
    battery_gwp_id = find_data(electric{1,1}.recharge_location,'gwp');
    battery_source_id = find_data(electric{1,1}.source_type,'source');
    vehicle.components{electric_id}.emissions.recharge = (electric{1}.mass*databases.source.specific_energy{battery_source_id} - electric{1}.remaining_energy)*databases.gwp.value{battery_gwp_id}*0.001/(3.6*10^6);
    vehicle.components{electric_id}.emissions.production = electric{1}.mass*databases.source.production_EI_CO2_per_kwh{battery_source_id}*databases.source.specific_energy{battery_source_id}/(3.6e6)/databases.source.max_num_cycles{battery_source_id};
end

if ~isempty(hydrogen) % emissions during production
    hydrogen_source_id = find_data(hydrogen{1,1}.source_type,'source');
    [fuel_cell,fuel_cell_id] = find_by_type(vehicle.components, 'fuel_cell');
    fuel_cell_source_id = find_data(fuel_cell{1,1}.component_type,'source');
    vehicle.components{hydrogen_id}.emissions.production = hydrogen{1}.mass*databases.source.specific_energy{hydrogen_source_id}*databases.source.production_EI_CO2_per_kwh{hydrogen_source_id}/(3.6*10^6);
    vehicle.components{fuel_cell_id}.emissions.production = fuel_cell{1}.max_power*databases.source.production_EI_CO2_per_kw{fuel_cell_source_id}*mission.time/1000/3600/databases.source.max_hours_operation{fuel_cell_source_id};
end








