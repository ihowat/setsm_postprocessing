from .create_neighbor_index import create_neighbor_index
from .create_working_dirs import create_working_dirs
from .link_matfiles_to_stage import link_files_to_stage
from .link_source_matfiles import link_source_matfiles
from .mosaic_outputs import mosaic_outputs
from .show_settings import show_settings
from .wrap_batch_merge_tile_buffer import merge_buffers
from .wrap_batch_register_tiles import (
    coreg_matfiles,
    water_flatten_matfiles,
)
from .wrap_batch_tiles2tif_v4 import coreg_debug, export_final_tifs

COMMANDS = [
    show_settings,
    create_working_dirs,
    link_source_matfiles,
    coreg_debug,
    link_files_to_stage,
    mosaic_outputs,
    coreg_matfiles,
    water_flatten_matfiles,
    create_neighbor_index,
    merge_buffers,
    export_final_tifs,
]
