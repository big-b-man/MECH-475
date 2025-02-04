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

function rotor = balance_forces_vertical(rotor, vehicle_mass)
global constants
denominator = 0;
for i = 1 : length(rotor) 
    denominator = denominator + rotor{i}.number * rotor_area(rotor{i}) / rotor_area(rotor{1});
end

for i = 1 : length(rotor)
    if i == 1
        rotor{i}.thrust = vehicle_mass * constants.g / denominator;
    else
        rotor{i}.thrust = rotor_area(rotor{i})/rotor_area(rotor{1}) * rotor{1}.thrust;
    end
end