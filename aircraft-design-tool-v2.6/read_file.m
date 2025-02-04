function data = read_file(filename,type)

temp = readcell(filename);

if strcmp(type,'source')
    data.type = temp(2:size(temp,1),1);
    data.specific_energy = temp(2:size(temp,1),2);
    data.operation_EI_CO2 = temp(2:size(temp,1),3);
    data.production_EI_CO2_per_kwh = temp(2:size(temp,1),4);
    data.production_EI_CO2_per_kw = temp(2:size(temp,1),5);
    data.max_hours_operation = temp(2:size(temp,1),6);
    data.max_num_cycles = temp(2:size(temp,1),7);
end

if strcmp(type,'gwp')
    data.region = temp(2:size(temp,1),1);
    data.value = temp(2:size(temp,1),2);
end