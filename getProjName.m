function [projabbrev, projstr] = getProjName(tileName,projstr)

    projstr_arctic_proper = 'polar stereo north';
    projstr_antarctic_proper = 'polar stereo south';

    projstr_arctic_accepted = {
        projstr_arctic_proper,
        'psn',
        'arctic',
        'arcticdem',
        'epsg:3413',
        '3413'
    };

    projstr_antarctic_accepted = {
        projstr_antarctic_proper,
        'pss',
        'antarctic',
        'rema',
        'epsg:3031',
        '3031'
    };

    projabbrev = '';

    if isempty(tileName) && isempty(projstr)
        error("Both 'tileName' and 'projstr' arguments cannot be empty");

    elseif ~isempty(tileName) && startsWith(tileName,'utm')
        tileName_parts = strsplit(tileName,'_');
        tileProjName = tileName_parts{1};
        projabbrev = tileProjName;
        projstr = tileProjName;

    elseif isempty(projstr)
%        error("Projection is undefined for tileName '%s'", tileName);
        return;

    elseif any(strcmpi(projstr_arctic_accepted, projstr))
        projabbrev = 'psn';
        projstr = projstr_arctic_proper;

    elseif any(strcmpi(projstr_antarctic_accepted, projstr))
        projabbrev = 'pss';
        projstr = projstr_antarctic_proper;

    else
        error("Projection is undefined for tileName '%s' and projstr '%s'", tileName, projstr);
    end

end
