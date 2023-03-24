function writeDataFiles()
    files = {dir("../data/*.mat").name};
    [~, nfiles] = size(files);

    pleth = zeros(144001, nfiles);
    co2 = zeros(144001, nfiles);

    for i = 1:nfiles
        dt = load("../data/" + cell2mat(files(i)));
        pleth(:, i) = dt.signal.pleth.y;
        co2(:, i) = dt.signal.co2.y;
    end

    Tpleth = array2table(pleth);
    Tpleth.Properties.VariableNames(:) = files;

    Tco2 = array2table(co2);
    Tco2.Properties.VariableNames(:) = files;

    writetable(Tpleth,'../data/pleth_signals.csv');
    writetable(Tco2,'../data/co2_signals.csv');

end