
import math


class InvalidArgumentError(Exception):
    def __init__(self, msg=""):
        super(Exception, self).__init__(msg)
class UnsupportedMethodError(Exception):
    def __init__(self, msg=""):
        super(Exception, self).__init__(msg)
class TileNameError(Exception):
    def __init__(self, msg=""):
        super(Exception, self).__init__(msg)


def parse_supertile_name(tile_name):
    tile_parts = tile_name.split('_')
    mosaic, super_row, super_col = None, None, None

    tile_name_valid = True

    if len(tile_parts) == 3:
        utm_zone, super_row, super_col = tile_parts
        mosaic = utm_zone
    elif len(tile_parts) == 2:
        super_row, super_col = tile_parts
        mosaic = 'polar'
    else:
        tile_name_valid = False

    for num in [super_row, super_col]:
        if not num.isdigit():
            tile_name_valid = False
            break

    if not tile_name_valid:
        raise TileNameError(
            "Invalid supertile name: '{}'".format(tile_name)
        )

    return mosaic, super_row, super_col


def group_supertiles_all_together(tile_list):
    group_tiles_dict = {}
    group_args_dict = None

    for supertile in tile_list:
        mosaic, super_row, super_col = parse_supertile_name(supertile)
        group_key = '{}_XX_YY'.format(mosaic)
        if group_key not in group_tiles_dict:
            group_tiles_dict[group_key] = []
        group_tiles_dict[group_key].append(supertile)

    return group_tiles_dict, group_args_dict


def group_tiles_all_separate(tile_list):
    group_tiles_dict = {tile: [tile] for tile in tile_list}
    group_args_dict = None
    return group_tiles_dict, group_args_dict


def group_supertiles_by_block(tile_list, block_size=3):
    group_tiles_dict = {}
    group_args_dict = None

    for supertile in tile_list:

        mosaic, super_row, super_col = parse_supertile_name(supertile)

        super_row = int(super_row)
        super_col = int(super_col)

        row_key = '{:02d}'.format(int(math.floor(max(0, super_row-1) / float(block_size))*block_size + 1))
        col_key = '{:02d}'.format(int(math.floor(max(0, super_col-1) / float(block_size))*block_size + 1))

        group_key = '{}_{}_{}'.format(mosaic, row_key, col_key)

        if group_key not in group_tiles_dict:
            group_tiles_dict[group_key] = []
        group_tiles_dict[group_key].append(supertile)

    return group_tiles_dict, group_args_dict


def group_supertiles_by_stripe(
        tile_list,
        stripe_offset,
        stripe_gap=2,
        dimension='col'
):
    arg_dimension_choices = ['row', 'col']
    if dimension not in arg_dimension_choices:
        raise InvalidArgumentError(
            "`dimension` argument must be one of {}, but was {}".format(
                arg_dimension_choices, dimension
            )
        )

    group_tiles_dict = {}
    group_args_dict = None

    for supertile in tile_list:

        mosaic, super_row, super_col = parse_supertile_name(supertile)

        super_row = int(super_row)
        super_col = int(super_col)

        check_tile_idx, row_key, col_key = None, None, None
        if dimension == 'row':
            check_tile_idx = super_row
            row_key = '{:02d}'.format(check_tile_idx)
            col_key = 'XX'
        elif dimension == 'col':
            check_tile_idx = super_col
            row_key = 'XX'
            col_key = '{:02d}'.format(check_tile_idx)

        if (check_tile_idx - stripe_offset) % (stripe_gap+1) == 0:
            pass
        else:
            continue

        group_key = '{}_{}_{}'.format(mosaic, row_key, col_key)

        if group_key not in group_tiles_dict:
            group_tiles_dict[group_key] = []
        group_tiles_dict[group_key].append(supertile)

    return group_tiles_dict, group_args_dict


def group_supertiles_by_quad_stripes_2gap(
        tile_list,
        stripe_offset,
        dimension='col'
):
    arg_dimension_choices = ['row', 'col']
    if dimension not in arg_dimension_choices:
        raise InvalidArgumentError(
            "`dimension` argument must be one of {}, but was {}".format(
                arg_dimension_choices, dimension
            )
        )

    group_tiles_dict = {}
    group_args_dict = {}

    for supertile in tile_list:

        mosaic, super_row, super_col = parse_supertile_name(supertile)

        super_row = int(super_row)
        super_col = int(super_col)

        check_supertile_idx, super_row_key, super_col_key = None, None, None
        if dimension == 'row':
            check_supertile_idx = super_row
            super_row_key = '{:02d}'.format(check_supertile_idx)
            super_col_key = 'XX'
        elif dimension == 'col':
            check_supertile_idx = super_col
            super_row_key = 'XX'
            super_col_key = '{:02d}'.format(check_supertile_idx)

        if (check_supertile_idx + stripe_offset + 1) % 3 == 0:
            continue
        else:
            pass

        if (check_supertile_idx + stripe_offset) % 3 == 0:
            quad_idx = 1
        else:
            quad_idx = 2

        quad_row_key, quad_col_key, quad_args = None, None, None
        if dimension == 'row':
            quad_row_key = str(quad_idx)
            quad_col_key = 'Y'
            quad_args = "'quad_row',{}".format(quad_idx)
        elif dimension == 'col':
            quad_row_key = 'Y'
            quad_col_key = str(quad_idx)
            quad_args = "'quad_col',{}".format(quad_idx)

        group_key = '_'.join([
            mosaic,
            super_row_key, super_col_key,
            quad_row_key, quad_col_key
        ])

        if group_key not in group_tiles_dict:
            group_tiles_dict[group_key] = []
        group_tiles_dict[group_key].append(supertile)

        if group_key not in group_args_dict:
            group_args_dict[group_key] = quad_args
        else:
            assert group_args_dict[group_key] == quad_args

    return group_tiles_dict, group_args_dict


def group_supertiles_by_quad_stripes_0gap(
        tile_list,
        dimension='col'
):
    arg_dimension_choices = ['row', 'col']
    if dimension not in arg_dimension_choices:
        raise InvalidArgumentError(
            "`dimension` argument must be one of {}, but was {}".format(
                arg_dimension_choices, dimension
            )
        )

    group_tiles_dict = {}
    group_args_dict = {}

    for supertile in tile_list:

        mosaic, super_row, super_col = parse_supertile_name(supertile)

        super_row = int(super_row)
        super_col = int(super_col)

        check_supertile_idx, super_row_key, super_col_key = None, None, None
        if dimension == 'row':
            check_supertile_idx = super_row
            super_row_key = '{:02d}'.format(check_supertile_idx)
            super_col_key = 'XX'
        elif dimension == 'col':
            check_supertile_idx = super_col
            super_row_key = 'XX'
            super_col_key = '{:02d}'.format(check_supertile_idx)

        for quad_idx in (1, 2):

            quad_row_key, quad_col_key, quad_args = None, None, None
            if dimension == 'row':
                quad_row_key = str(quad_idx)
                quad_col_key = 'Y'
                quad_args = "'quad_row',{}".format(quad_idx)
            elif dimension == 'col':
                quad_row_key = 'Y'
                quad_col_key = str(quad_idx)
                quad_args = "'quad_col',{}".format(quad_idx)

            group_key = '_'.join([
                mosaic,
                super_row_key, super_col_key,
                quad_row_key, quad_col_key
            ])

            if group_key not in group_tiles_dict:
                group_tiles_dict[group_key] = []
            group_tiles_dict[group_key].append(supertile)

            if group_key not in group_args_dict:
                group_args_dict[group_key] = quad_args
            else:
                assert group_args_dict[group_key] == quad_args

    return group_tiles_dict, group_args_dict


def group_supertiles_by_row_or_col(
        tile_list,
        dimension,
        resolution,
):
    arg_dimension_choices = ['row', 'col']
    if dimension not in arg_dimension_choices:
        raise InvalidArgumentError(
            "`dimension` argument must be one of {}, but was {}".format(
                arg_dimension_choices, dimension
            )
        )

    arg_resolution_choices = [10, 2]
    if resolution not in arg_resolution_choices:
        raise InvalidArgumentError(
            "`resolution` argument must be one of {}, but was {}".format(
                arg_resolution_choices, resolution
            )
        )

    if resolution == 10:
        return group_supertiles_by_stripe(tile_list, 0, 0, dimension)
    elif resolution == 2:
        return group_supertiles_by_quad_stripes_0gap(tile_list, dimension)


def group_supertiles_for_processing(script_args, tile_list):
    group_tiles_dict, group_args_dict = None, None

    if not hasattr(script_args, 'process_group') or script_args.process_group is None:
        process_group = 'all'
    else:
        process_group = script_args.process_group

    resolution = script_args.resolution

    if process_group == 'all':
        group_tiles_dict, group_args_dict = group_supertiles_all_together(tile_list)

    elif process_group == 'separate':
        group_tiles_dict, group_args_dict = group_tiles_all_separate(tile_list)

    elif process_group == 'block':
        group_tiles_dict, group_args_dict = group_supertiles_by_block(tile_list)

    elif process_group.startswith('stripe-'):
        stripe_offset = int(process_group.split('-')[1]) - 1
        if resolution == 2:
            group_tiles_dict, group_args_dict = group_supertiles_by_quad_stripes_2gap(tile_list, stripe_offset)
        elif resolution == 10:
            group_tiles_dict, group_args_dict = group_supertiles_by_stripe(tile_list, stripe_offset)

    elif process_group in ('row', 'column'):
        if process_group == 'row':
            dimension = 'row'
        elif process_group == 'column':
            dimension = 'col'
        group_tiles_dict, group_args_dict = group_supertiles_by_row_or_col(tile_list, dimension, resolution)

    else:
        raise UnsupportedMethodError(
            "Script argument --process-group='{}' method is not supported".format(process_group)
        )

    return group_tiles_dict, group_args_dict
