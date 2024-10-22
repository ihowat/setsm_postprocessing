from . import (
    build_subtiles,
    merge_buffers,
    add_land_mask,
    export_products,
    show_settings,
)

COMMANDS = [
    build_subtiles.init_bst,
    build_subtiles.run_bst,
    merge_buffers.create_neighbor_index,
    merge_buffers.backup_matfiles,
    merge_buffers.merge_tile_buffers,
    add_land_mask.add_land_mask,
    export_products.export_products,
    export_products.build_vrts,
    show_settings.show_settings,
]
