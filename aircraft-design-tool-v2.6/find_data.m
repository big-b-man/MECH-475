function id = find_data(name,database)
global databases

if strcmp(database,'source')
    for i = 1:size(databases.source.type,1)
        if strcmp(databases.source.type{i,1},name)
            id = i;
        end
    end
elseif strcmp(database,'gwp')
    for i = 1:size(databases.gwp.region,1)
        if strcmp(databases.gwp.region{i,1},name)
            id = i;
        end
    end
end