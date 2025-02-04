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

function [elems, ids] = find_network_components(vehicle, network)
elems = {};
ids = [];
for i = 1 : length(network.layout)
    for k = 1 : length(network.layout{i})
        for j = 1 : length(vehicle.components)
            if strcmp(network.layout{i}(k).name, vehicle.components{j}.name)
                e = vehicle.components{j};
                fields = fieldnames(network.layout{i}(k));
                for m = 1 : length(fields)
                    e.(fields{m}) = network.layout{i}(k).(fields{m});
                end
                elems{i,k} = e;
                ids(i,k) = j;
                break;
            end
        end
    end
end
