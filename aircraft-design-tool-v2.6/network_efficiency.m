
% Compute network efficiencies (in series)
function e = network_efficiency(network, type)

id = 0;

% Determine location of component to start calculating efficiency
for i = 1 : length(network)
    for j = 1 : size(network,2)
        if (~isempty(network{i,j})) && is_type(network{i,j}, type)
            id = i;
            break
        end
    end
end

% Initialize matrix with efficiencies
e = zeros(size(network,2),1);
e(:,1) = 1.0;

for i = 1 : length(network)
    if i > id
        for j = 1 : size(network,2)
            if (~isempty(network{i,j})) && isfield(network{i,j},'efficiency')  
                e(j,1) = e(j,1) * network{i,j}.efficiency;
            elseif isempty(network{i,j}) 
                for k = (j - 1) : 1
                    if (~isempty(network{i,k})) && isfield(network{i,k},'efficiency')
                        e(j,1) = e(j,1) * network{i,k}.efficiency;
                    end
                end
            end
        end
    end
end
