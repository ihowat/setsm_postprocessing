function [parsed_tilename, supername, quadname] = parseTileName(tilename)

allow_suffix = false;

prefix = [];
supertile = [];
quadtile = [];
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
end
if isfield(tokens, 'quadtile') && ~isempty(tokens.quadtile)
    quadtile = strip(tokens.quadtile, 'left', '_');
end
if isfield(tokens, 'suffix') && ~isempty(tokens.suffix)
    suffix = tokens.suffix;
end

tilename_parts = {prefix, supertile, quadtile};
tilename_parts(cellfun(@(s) isempty(s), tilename_parts)) = [];
parsed_tilename = strjoin(tilename_parts, '_');

supername_parts = {prefix, supertile};
supername_parts(cellfun(@(s) isempty(s), supername_parts)) = [];
supername = strjoin(supername_parts, '_');

quadname_parts = {prefix, supertile, quadtile};
quadname_parts(cellfun(@(s) isempty(s), quadname_parts)) = [];
quadname = strjoin(quadname_parts, '_');
