function [datestring_yyyymmdd] = parsePairnameDatestring(string)

prefix = [];
pairname = [];
sensor = [];
datestring_yyyymmdd = [];
catalogid1 = [];
catalogid2 = [];
suffix = [];

pairname_re = '^(?<prefix>.*)?(?<sensor>[A-Z][A-Z0-9]{2}[0-9])_(?<date>[0-9]{8})_(?<catid1>[0-9A-F]{16})_(?<catid2>[0-9A-F]{16})(?<suffix>.*)?$';

[tokens, match_idx] = regexp(string, pairname_re, 'names');
if isempty(match_idx)
    error("Cannot parse pairname parts with regex '%s' from input tilename string: %s", pairname_re, string);
end

%if isfield(tokens, 'prefix')
%    prefix = tokens.prefix;
%end
%sensor = tokens.sensor;
datestring_yyyymmdd = tokens.date;
%catalogid1 = tokens.catid1;
%catalogid2 = tokens.catid2;
%if isfield(tokens, 'suffix')
%    suffix = tokens.suffix;
%end
%
%pairname_parts = [sensor, datestring_yyyymmdd, catalogid1, catalogid2]
%pairname = strjoin(pairname_parts, '_');
