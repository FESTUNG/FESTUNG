pd.configSource = 'rotation';
pd.gridSource = 'square';

%% 1 cylinder
r = 0.0225 / 4;
for radius = 1:4
    pd.transportData.cH0Cont{1} = @(x1, x2) ( ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= r) ) .* 0.002;
    pd.transportData.isMask = true;
    main('swe_transport',pd);
    pd.transportData.isMask = false;
    main('swe_transport',pd);
    r = r * sqrt(2);
end
    

%% 4 cylinders
r = 0.0225 / 4;
for radius = 1:4
    pd.transportData.cH0Cont{1} = @(x1, x2)( ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= r) ...
                                           + ((x1 - 0.25).^2 + (x2 - 0.5).^2 <= r) ...
                                           + ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= r) ...
                                           + ((x1 - 0.75).^2 + (x2 - 0.5).^2 <= r) ) .* 0.002;
    pd.transportData.isMask = true;
    main('swe_transport',pd);
    pd.transportData.isMask = false;
    main('swe_transport',pd);
    r = r * sqrt(2);
end