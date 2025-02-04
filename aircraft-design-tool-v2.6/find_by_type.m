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

function [elems, ids] = find_by_type(data, type)
elems = {};
ids = [];
count = 0;
type_tags = split(type, '.');
for i = 1 : length(data)
    for k = 1 : size(data,2)
        if (~isempty(data{i,k}))
            match = true;

            elem_tags = split(data{i,k}.type, '.');
            for j = 1 : length(type_tags)
                if ~strcmp(elem_tags{j}, type_tags{j})
                    match = false;
                    break;
                end
            end
    
            if (match == true)
                count = count + 1;
                elems{count,1} = data{i,k};
                ids(count,1) = i;
            end
        end
    end
end