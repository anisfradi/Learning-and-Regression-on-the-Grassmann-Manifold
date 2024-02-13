function var = computeVar( omas, omas_basis, flag )
% flag: 1 -- Frechet mean, 2 -- within sample mean, 3 -- direct diff 

switch flag
    case 1
        fMean = grfmean( omas, 1e-6 );
        sum = 0;
        for iI = 1:length(omas_basis)
            sum = sum + grarc( fMean, omas_basis{iI} );
        end
        var = sum / length(omas_basis);
    case 2
        minValue = 100000;
        for iI = 1:length(omas)
            sum = 0;
            for iJ = 1:length(omas_basis)
                sum = sum + grarc(omas{iI}, omas_basis{iJ});
            end
            sum = sum / length(omas_basis);
            if iI == 1 || sum < minValue
                minValue = sum;
            end
        end
        var = minValue;
    case 3
        sum = 0;
        for iI = 1:length(omas)
            sum = sum + grarc(omas{iI}, omas_basis{iI});
        end
        var = sum/length(omas);
    otherwise
        error( 'Unknown input' );
end