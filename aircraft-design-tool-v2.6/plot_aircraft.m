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

function vehicle = plot_aircraft(vehicle)

figure;
legend;
hold on;
a = gca;
a.XLabel.String = 'x';
a.YLabel.String = "y";

weighed_cg = 0;

% Plot fuselage
[fuselage, fuselage_id] = find_by_type(vehicle.components, 'fuselage');

% Fuselage
for i = 1:length(fuselage{1,1}.diameters)
    pos(i,1) = fuselage{1,1}.length/(length(fuselage{1,1}.diameters)-1)*(i-1);
    pos(i,2) = -fuselage{1,1}.diameters(i)/2;
end

sum_diameters = 0;
sum_cg_pos = 0;
for i = 1:length(fuselage{1,1}.diameters)
    sum_diameters = sum_diameters + fuselage{1,1}.diameters(i);
    sum_cg_pos = sum_cg_pos + pos(i,1) * fuselage{1,1}.diameters(i);
end

vehicle.components{fuselage_id}.position.cg = sum_cg_pos/sum_diameters;
weighed_cg = weighed_cg + vehicle.components{fuselage_id}.position.cg * fuselage{1,1}.mass;

for i = 1:length(fuselage{1,1}.diameters)
    pos(i+length(fuselage{1,1}.diameters),1) = fuselage{1,1}.length/(length(fuselage{1,1}.diameters)-1)*(length(fuselage{1,1}.diameters)-i);
    pos(i+length(fuselage{1,1}.diameters),2) = fuselage{1,1}.diameters(length(fuselage{1,1}.diameters)-(i-1))/2;
end

for i = 1:(length(pos)-1)
    plot([pos(i,1),pos(i+1,1)],[pos(i,2),pos(i+1,2)], "Color", 'b', 'HandleVisibility','off');
end

plot([pos(length(pos),1),pos(1,1)],[pos(length(pos),2),pos(1,2)], "Color", 'b', 'HandleVisibility','off');
scatter(vehicle.components{fuselage_id}.position.cg,0, 'b', 'DisplayName', 'C.G. Fuselage', 'Marker', '+')

% Booms
[boom, boom_id] = find_by_type(vehicle.components, 'boom');

if ~isempty(boom)
    for i = 1:size(vehicle.components{boom_id}.position.structure,1)
        vehicle.components{boom_id}.position.cg(i,:) = [(vehicle.components{boom_id}.position.structure(i,1,1) + vehicle.components{boom_id}.position.structure(i,2,1))/2,(vehicle.components{boom_id}.position.structure(i,1,2) + vehicle.components{boom_id}.position.structure(i,2,2))/2];
        weighed_cg = weighed_cg + vehicle.components{boom_id}.position.cg(i,1) * boom{1,1}.mass;
        plot([vehicle.components{boom_id}.position.structure(i,1,1),vehicle.components{boom_id}.position.structure(i,2,1)],[vehicle.components{boom_id}.position.structure(i,1,2),vehicle.components{boom_id}.position.structure(i,2,2)], "Color", 'b', 'HandleVisibility','off');
        scatter(vehicle.components{boom_id}.position.cg(i,1),vehicle.components{boom_id}.position.cg(i,2), 'r', 'DisplayName', 'C.G. Boom', 'Marker', 'v')
    end
end

if ~is_type(find_by_type(vehicle.components,'aircraft'),'aircraft.rotary_wing')
    % Main wing
    [wing, ~] = find_by_type(vehicle.components, 'wing.main');
    wing{1,1}.position.cg = wing{1,1}.position.leading_edge + (wing{1,1}.span/6)*((1+2*wing{1,1}.taper_ratio)/(1+wing{1,1}.taper_ratio))*tand(wing{1,1}.sweep_le)+wing{1,1}.airfoil.xc_max*wing{1,1}.mean_chord;
    wing{1,1}.position.ac = (wing{1,1}.span/6)*((1+2*wing{1,1}.taper_ratio)/(1+wing{1,1}.taper_ratio))*tand(wing{1,1}.sweep_le) + wing{1,1}.mean_chord*0.25;
    weighed_cg = weighed_cg + wing{1,1}.position.cg * wing{1,1}.mass;

    p1 = [wing{1,1}.position.leading_edge,0];
    p2 = [wing{1,1}.position.leading_edge + (wing{1,1}.span/2 * tand(wing{1,1}.sweep_le)), -wing{1,1}.span/2];
    p3 = [wing{1,1}.position.leading_edge + (wing{1,1}.span/2 * tand(wing{1,1}.sweep_le)) + wing{1,1}.taper_ratio * wing{1,1}.root_chord, -wing{1,1}.span/2];
    p4 = [wing{1,1}.position.leading_edge + wing{1,1}.root_chord,0];
    p5 = [wing{1,1}.position.leading_edge + (wing{1,1}.span/2 * tand(wing{1,1}.sweep_le)) + wing{1,1}.taper_ratio * wing{1,1}.root_chord, wing{1,1}.span/2];
    p6 = [wing{1,1}.position.leading_edge + (wing{1,1}.span/2 * tand(wing{1,1}.sweep_le)), wing{1,1}.span/2];
    p7 = [wing{1,1}.position.leading_edge,0];

    plot([p1(1,1),p2(1,1)],[p1(1,2),p2(1,2)], "Color", 'b', 'HandleVisibility','off');
    plot([p2(1,1),p3(1,1)],[p2(1,2),p3(1,2)], "Color", 'b', 'HandleVisibility','off');
    plot([p3(1,1),p4(1,1)],[p3(1,2),p4(1,2)], "Color", 'b', 'HandleVisibility','off');
    plot([p4(1,1),p5(1,1)],[p4(1,2),p5(1,2)], "Color", 'b', 'HandleVisibility','off');
    plot([p5(1,1),p6(1,1)],[p5(1,2),p6(1,2)], "Color", 'b', 'HandleVisibility','off');
    plot([p6(1,1),p7(1,1)],[p6(1,2),p7(1,2)], "Color", 'b', 'HandleVisibility','off');
    scatter(wing{1,1}.position.cg,0, 'r', 'DisplayName', 'C.G. Wing', 'Marker', '*');

    % Secondary wing (tandem)
    [sec_wing, ~] = find_by_type(vehicle.components, 'wing.secondary');
    if ~isempty(sec_wing)
        sec_wing{1,1}.position.cg = sec_wing{1,1}.position.leading_edge + (sec_wing{1,1}.span/6)*((1+2*sec_wing{1,1}.taper_ratio)/(1+sec_wing{1,1}.taper_ratio))*tand(sec_wing{1,1}.sweep_le)+sec_wing{1,1}.airfoil.xc_max*sec_wing{1,1}.mean_chord;
        sec_wing{1,1}.position.ac = sec_wing{1,1}.position.leading_edge + (sec_wing{1,1}.span/6)*((1+2*sec_wing{1,1}.taper_ratio)/(1+sec_wing{1,1}.taper_ratio))*tand(sec_wing{1,1}.sweep_le) + sec_wing{1,1}.mean_chord*0.25;
        weighed_cg = weighed_cg + sec_wing{1,1}.position.cg * sec_wing{1,1}.mass;

        p1 = [sec_wing{1,1}.position.leading_edge,0];
        p2 = [sec_wing{1,1}.position.leading_edge + (sec_wing{1,1}.span/2 * tand(sec_wing{1,1}.sweep_le)), -sec_wing{1,1}.span/2];
        p3 = [sec_wing{1,1}.position.leading_edge + (sec_wing{1,1}.span/2 * tand(sec_wing{1,1}.sweep_le)) + sec_wing{1,1}.taper_ratio * sec_wing{1,1}.root_chord, -sec_wing{1,1}.span/2];
        p4 = [sec_wing{1,1}.position.leading_edge + sec_wing{1,1}.root_chord,0];
        p5 = [sec_wing{1,1}.position.leading_edge + (sec_wing{1,1}.span/2 * tand(sec_wing{1,1}.sweep_le)) + sec_wing{1,1}.taper_ratio * sec_wing{1,1}.root_chord, sec_wing{1,1}.span/2];
        p6 = [sec_wing{1,1}.position.leading_edge + (sec_wing{1,1}.span/2 * tand(sec_wing{1,1}.sweep_le)), sec_wing{1,1}.span/2];
        p7 = [sec_wing{1,1}.position.leading_edge,0];

        plot([p1(1,1),p2(1,1)],[p1(1,2),p2(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p2(1,1),p3(1,1)],[p2(1,2),p3(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p3(1,1),p4(1,1)],[p3(1,2),p4(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p4(1,1),p5(1,1)],[p4(1,2),p5(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p5(1,1),p6(1,1)],[p5(1,2),p6(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p6(1,1),p7(1,1)],[p6(1,2),p7(1,2)], "Color", 'b', 'HandleVisibility','off');
        scatter(sec_wing{1,1}.position.cg,0, 'r', 'DisplayName', 'C.G. Secondary Wing', 'Marker', '*');
    end

    % Canard
    [canard, canard_id] = find_by_type(vehicle.components, 'wing.canard');
    if ~isempty(canard)
        vehicle.components{canard_id}.position.cg = canard{1,1}.position.leading_edge + (canard{1,1}.span/6)*((1+2*canard{1,1}.taper_ratio)/(1+canard{1,1}.taper_ratio))*tand(canard{1,1}.sweep_le)+canard{1,1}.airfoil.xc_max*canard{1,1}.mean_chord;
        weighed_cg = weighed_cg + vehicle.components{canard_id}.position.cg * canard{1,1}.mass;

        p1 = [canard{1,1}.position.leading_edge,0];
        p2 = [canard{1,1}.position.leading_edge + (canard{1,1}.span/2 * tand(canard{1,1}.sweep_le)), -canard{1,1}.span/2];
        p3 = [canard{1,1}.position.leading_edge + (canard{1,1}.span/2 * tand(canard{1,1}.sweep_le)) + canard{1,1}.taper_ratio * canard{1,1}.root_chord, -canard{1,1}.span/2];
        p4 = [canard{1,1}.position.leading_edge + canard{1,1}.root_chord,0];
        p5 = [canard{1,1}.position.leading_edge + (canard{1,1}.span/2 * tand(canard{1,1}.sweep_le)) + canard{1,1}.taper_ratio * canard{1,1}.root_chord, canard{1,1}.span/2];
        p6 = [canard{1,1}.position.leading_edge + (canard{1,1}.span/2 * tand(canard{1,1}.sweep_le)), canard{1,1}.span/2];
        p7 = [canard{1,1}.position.leading_edge,0];
    
        plot([p1(1,1),p2(1,1)],[p1(1,2),p2(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p2(1,1),p3(1,1)],[p2(1,2),p3(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p3(1,1),p4(1,1)],[p3(1,2),p4(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p4(1,1),p5(1,1)],[p4(1,2),p5(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p5(1,1),p6(1,1)],[p5(1,2),p6(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p6(1,1),p7(1,1)],[p6(1,2),p7(1,2)], "Color", 'b', 'HandleVisibility','off');
        scatter(vehicle.components{canard_id}.position.cg,0, 'r', 'DisplayName', 'C.G. Canard', 'Marker', '^');
    end

    % Horizontal Tail
    [htail, htail_id] = find_by_type(vehicle.components, 'wing.htail');
    if ~isempty(htail)
        vehicle.components{htail_id}.position.cg = htail{1,1}.position.leading_edge + (htail{1,1}.span/6)*((1+2*htail{1,1}.taper_ratio)/(1+htail{1,1}.taper_ratio))*tand(htail{1,1}.sweep_le)+htail{1,1}.airfoil.xc_max*htail{1,1}.mean_chord;
        weighed_cg = weighed_cg + vehicle.components{htail_id}.position.cg * htail{1,1}.mass;

        p1 = [htail{1,1}.position.leading_edge,0];
        p2 = [htail{1,1}.position.leading_edge + (htail{1,1}.span/2 * tand(htail{1,1}.sweep_le)), -htail{1,1}.span/2];
        p3 = [htail{1,1}.position.leading_edge + (htail{1,1}.span/2 * tand(htail{1,1}.sweep_le)) + htail{1,1}.taper_ratio * htail{1,1}.root_chord, -htail{1,1}.span/2];
        p4 = [htail{1,1}.position.leading_edge + htail{1,1}.root_chord,0];
        p5 = [htail{1,1}.position.leading_edge + (htail{1,1}.span/2 * tand(htail{1,1}.sweep_le)) + htail{1,1}.taper_ratio * htail{1,1}.root_chord, htail{1,1}.span/2];
        p6 = [htail{1,1}.position.leading_edge + (htail{1,1}.span/2 * tand(htail{1,1}.sweep_le)), htail{1,1}.span/2];
        p7 = [htail{1,1}.position.leading_edge,0];

        plot([p1(1,1),p2(1,1)],[p1(1,2),p2(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p2(1,1),p3(1,1)],[p2(1,2),p3(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p3(1,1),p4(1,1)],[p3(1,2),p4(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p4(1,1),p5(1,1)],[p4(1,2),p5(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p5(1,1),p6(1,1)],[p5(1,2),p6(1,2)], "Color", 'b', 'HandleVisibility','off');
        plot([p6(1,1),p7(1,1)],[p6(1,2),p7(1,2)], "Color", 'b', 'HandleVisibility','off');
        scatter(vehicle.components{htail_id}.position.cg,0, 'r', 'DisplayName', 'C.G. Horizontal Tail', 'Marker', '^');
    end

    % V-tail
    [vtail, vtail_id] = find_by_type(vehicle.components, 'wing.vtail');
    if ~isempty(vtail)
        vehicle.components{vtail_id}.position.cg = vtail{1,1}.position.leading_edge + (vtail{1,1}.half_span*2/6)*((1+2*vtail{1,1}.taper_ratio)/(1+vtail{1,1}.taper_ratio))*tand(vtail{1,1}.sweep_le)+vtail{1,1}.airfoil.xc_max*vtail{1,1}.mean_chord;
        weighed_cg = weighed_cg + vehicle.components{vtail_id}.position.cg * vtail{1,1}.mass;
        scatter(vehicle.components{vtail_id}.position.cg,0, 'r', 'DisplayName', 'C.G. Vertical Tail', 'Marker', 'x');
    end
end

% Other components
for i = 1 : length(vehicle.components)
    if is_type(vehicle.components{i}, 'fuselage') || is_type(vehicle.components{i}, 'wing') || is_type(vehicle.components{i}, 'boom')
        % Nothing
    elseif is_type(vehicle.components{i}, 'driver')
        marker = 'diamond';
        if isfield(vehicle.components{i}.position.cg, 'vtol')
            for j = 1:size(vehicle.components{i}.position.cg.vtol,1)
                weighed_cg = weighed_cg + vehicle.components{i}.position.cg.vtol(j,1) * vehicle.components{i}.mass;
                scatter(vehicle.components{i}.position.cg.vtol(j,1),vehicle.components{i}.position.cg.vtol(j,2), 'filled', 'DisplayName', vehicle.components{i}.name, 'Marker', marker);
                viscircles([vehicle.components{i}.position.cg.vtol(j,1),vehicle.components{i}.position.cg.vtol(j,2)],vehicle.components{i}.radius,'Color','b'); % draw circled rotors
            end
        elseif isfield(vehicle.components{i}.position.cg, 'forward')
            for j = 1:size(vehicle.components{i}.position.cg.forward,1)
                weighed_cg = weighed_cg + vehicle.components{i}.position.cg.forward(j,1) * vehicle.components{i}.mass;
                scatter(vehicle.components{i}.position.cg.forward(j,1),vehicle.components{i}.position.cg.forward(j,2), 'filled', 'DisplayName', vehicle.components{i}.name, 'Marker', marker);
            end
        end
    else
        if (~is_type(vehicle.components{i}, 'mass.empty'))
            for j = 1:size(vehicle.components{i}.position.cg,1)
                weighed_cg = weighed_cg + vehicle.components{i}.position.cg(j,1) * vehicle.components{i}.mass;
                if is_type(vehicle.components{i}, 'energy')
                    marker = 'square';
                elseif is_type(vehicle.components{i}, 'engine') || is_type(vehicle.components{i}, 'motor')
                    marker = 'v';
                else 
                    marker = 'o';
                end
                scatter(vehicle.components{i}.position.cg(j,1),vehicle.components{i}.position.cg(j,2), 'filled', 'DisplayName', vehicle.components{i}.name, 'Marker', marker);
            end
        end
    end
end
vehicle.x_cg = weighed_cg/vehicle.mass;
scatter(vehicle.x_cg, 0 , 'filled', 'DisplayName', 'Global C.G.', 'Marker', 'pentagram');
axis equal










