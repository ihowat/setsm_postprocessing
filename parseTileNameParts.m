function [prefix, supertile, quadtile, suffix] = parseTileNameParts(tilename, allow_suffix)

prefix = [];
supertile = [];
quadtile = [];
suffix = [];

tilename_re = '(?<prefix>utm\d{2}[ns]_)?(?<supertile>\d{2}_\d{2}s?)(?<quadtile>_\d_\d)?';
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
end
if isfield(tokens, 'quadtile') && ~isempty(tokens.quadtile)
    quadtile = strip(tokens.quadtile, 'left', '_');
end
if isfield(tokens, 'suffix') && ~isempty(tokens.suffix)
    suffix = tokens.suffix;
end
