function [prefix, super_row, super_col, quad_row, quad_col, suffix] = parseTileNameAll(tilename, allow_suffix)

prefix = [];
super_row = [];
super_col = [];
quad_row = [];
quad_col = [];
suffix = [];

tilename_re = '(?<prefix>utm\d{2}[ns]_)?(?<supertile>\d{2}_\d{2})(?<quadtile>_\d_\d)?';
if allow_suffix
    tilename_re = ['^', tilename_re, '(?<suffix>.*$)'];
else
    tilename_re = ['^', tilename_re, '$'];
end

[tokens, match_idx] = regexp(tilename, tilename_re, 'names');
if isempty(match_idx)
    error("Cannot parse tilename parts with regex '%s' from input tilename string: %s", tilename_re, tilename);
end

if isfield(tokens, 'prefix') && ~isempty(tokens.prefix)
    prefix = strip(tokens.prefix, 'right', '_');
end
if isfield(tokens, 'supertile') && ~isempty(tokens.supertile)
    supertile = tokens.supertile;
    parts = split(supertile, '_');
    super_row = str2double(parts{1});
    super_col = str2double(parts{2});
end
if isfield(tokens, 'quadtile') && ~isempty(tokens.quadtile)
    quadtile = strip(tokens.quadtile, 'left', '_');
    parts = split(quadtile, '_');
    quad_row = str2double(parts{1});
    quad_col = str2double(parts{2});
end
if isfield(tokens, 'suffix') && ~isempty(tokens.suffix)
    suffix = tokens.suffix;
end
