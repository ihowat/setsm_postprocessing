function [epsg,utm_hemi,utm_zone] = getProjstrInfo(projstr)

epsg = [];
utm_hemi = [];
utm_zone = [];

if contains(projstr,'UTM','IgnoreCase',true) ...
    || ~isempty(regexp(projstr, '^\s*(?:UTM)?\s*\d{1,2}\s*(?:North|South|N|S)\s*$', 'ignorecase')) ...
    || ~isempty(regexp(projstr, '^\s*(?:UTM)?\s*(?:North|South|N|S)\s*\d{1,2}\s*$', 'ignorecase'))

    projstr_trim = regexprep(projstr,'utm','','ignorecase');

    if contains(projstr,'North','IgnoreCase',true)
        utm_hemi = 'north';
        projstr_trim = regexprep(projstr_trim,'north','','ignorecase');
    elseif contains(projstr,'South','IgnoreCase',true)
        utm_hemi = 'south';
        projstr_trim = regexprep(projstr_trim,'south','','ignorecase');
    elseif contains(projstr,'N','IgnoreCase',true)
        utm_hemi = 'north';
        projstr_trim = regexprep(projstr_trim,'n','','ignorecase');
    elseif contains(projstr,'S','IgnoreCase',true)
        utm_hemi = 'south';
        projstr_trim = regexprep(projstr_trim,'s','','ignorecase');
    else
        error('Cannot parse hemisphere information from UTM ''projstr'': %s', projstr);
    end

    utm_zone = str2num(projstr_trim);
    if isempty(utm_zone)
        error('Cannot parse zone number from UTM ''projstr'': %s', projstr);
    end

    if ~isempty(utm_hemi) && ~isempty(utm_zone)
        if strcmp(utm_hemi, 'north')
            epsg = 32600 + utm_zone;
        elseif strcmp(utm_hemi, 'south')
            epsg = 32700 + utm_zone;
        end
    end
    if isempty(epsg)
        error('Failed to determine EPSG code for UTM ''projstr'': %s', projstr);
    end

else
    switch lower(projstr)
        case 'polar stereo north'; epsg = 3413;
        case 'polar stereo south'; epsg = 3031;
    end
    if isempty(epsg)
        error('Failed to lookup EPSG code for ''projstr'': %s', projstr);
    end
end

if isempty(epsg)
    error('Failed to determine EPSG code for ''projstr'': %s', projstr);
end
