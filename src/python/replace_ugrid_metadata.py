#!/usr/bin/env python3
"""
Python utility that replaces UGRID metadata of an LFRic netCDF file with metadata from an
LFRic mesh file, and writes all metadata and field data to a new UGRID netCDF output file.
"""
import argparse as ap
import logging
import netCDF4 as nc


def get_fields(nc_file):
    """
    Scan input netCDF file and identify all variables that are associated with model data fields
    """

    fields = []
    for var in nc_file.variables.values():
        if all(att in dir(var) for att in ['mesh', 'units', 'location']):
            fields.append(var)

    logging.info('Found %i field variables.', len(fields))

    return fields


def get_vertical_axes(nc_file):
    """
    Scan input netCDF file and identify vertical axes
    """

    vertical_axes = []
    for var_name, var in nc_file.variables.items():
        if var_name in ('full_levels', 'half_levels'):
            vertical_axes.append(var)

    logging.info('Found %i vertical axes.', len(vertical_axes))

    return vertical_axes


def get_time_vars(nc_file):
    """
    Scan input netCDF file and identify time variables
    """

    time_vars = []
    for var_name, var in nc_file.variables.items():
        if var_name in ('time_instant', 'time_instant_bounds'):
            time_vars.append(var)

    logging.info('Found %i time variables.', len(time_vars))

    return time_vars


def get_ugrid_mesh_description(nc_file):
    """
    Collect UGRID mesh description from mesh file
    """

    # Look for UGRID mesh description dummy variable
    mesh_names = []
    for var_name, var in nc_file.variables.items():
        if 'cf_role' in dir(var):
            if var.cf_role == 'mesh_topology':
                mesh_names.append(var_name)

    if len(mesh_names) != 1:
        logging.critical('Expected 1 UGRID mesh description but found %i', len(mesh_names))
        raise RuntimeError()

    logging.info('Found UGRID mesh "%s"', mesh_names[0])

    ugrid_topology_var = nc_file.variables[mesh_names[0]]

    # Collect expected UGRID location and connectivity variables
    ugrid_varnames = set()
    ugrid_varnames.update(ugrid_topology_var.node_coordinates.split(" "))
    ugrid_varnames.update(ugrid_topology_var.face_coordinates.split(" "))
    ugrid_varnames.add(ugrid_topology_var.face_node_connectivity)
    ugrid_varnames.add(ugrid_topology_var.edge_node_connectivity)
    ugrid_varnames.add(ugrid_topology_var.face_edge_connectivity)
    ugrid_varnames.add(ugrid_topology_var.face_face_connectivity)

    ugrid_vars = [ugrid_topology_var]
    for varname in ugrid_varnames:
        ugrid_vars.append(nc_file.variables[varname])

    return ugrid_vars


def write_file(nc_file, ugrid_vars, vert_axes, time_vars, fields):
    """
    Write netCDF output file
    """

    #
    # UGRID
    #

    # Collect UGRID dimensions
    ugrid_dims = set()
    for ugrid_var in ugrid_vars:
        ugrid_dims.update(ugrid_var.get_dims())

    mesh_name = ugrid_vars[0].name
    output_mesh_name = 'Mesh2d'

    # Generate dictionary for translating dimension names
    ugrid_dim_translate = dict()
    for dim in ugrid_dims:
        ugrid_dim_translate[dim.name] = dim.name.replace(mesh_name, output_mesh_name)

    # Create UGRID dimensions in output file
    for dim in ugrid_dims:
        nc_file.createDimension(ugrid_dim_translate[dim.name], dim.size)

    # Create new field variables in output file
    for ugrid_var in ugrid_vars:
        output_var_dims = []
        for dim in ugrid_var.dimensions:
            output_var_dims.append(ugrid_dim_translate[dim])
        output_var = nc_file.createVariable(ugrid_var.name.replace(mesh_name, output_mesh_name),
                                            ugrid_var.datatype, output_var_dims)
        output_var[:] = ugrid_var[:]
        for att in ugrid_var.ncattrs():
            att_value = ugrid_var.getncattr(att)
            if isinstance(att_value, str):
                output_var.setncattr(att, att_value.replace(mesh_name, output_mesh_name))
            else:
                output_var.setncattr(att, att_value)

    #
    # Vertical axes
    #

    vert_dims = set()
    for vert_axis in vert_axes:
        vert_dims.update(vert_axis.get_dims())

    for dim in vert_dims:
        nc_file.createDimension(dim.name, dim.size)

    for vert_axis in vert_axes:
        output_axis = nc_file.createVariable(vert_axis.name, vert_axis.datatype,
                                             vert_axis.dimensions)
        output_axis[:] = vert_axis[:]
        for att in vert_axis.ncattrs():
            output_axis.setncattr(att, vert_axis.getncattr(att))

    #
    # Time axis
    #

    time_dims = set()
    for time_var in time_vars:
        time_dims.update(time_var.get_dims())

    for dim in time_dims:
        nc_file.createDimension(dim.name, dim.size)

    for time_var in time_vars:
        output_var = nc_file.createVariable(time_var.name, time_var.datatype, time_var.dimensions)
        output_var[:] = time_var[:]
        for att in time_var.ncattrs():
            output_var.setncattr(att, time_var.getncattr(att))

    #
    # Fields
    #

    # Collect dimensions from input file
    field_dims = set()
    for field in fields:
        field_dims.update(field.get_dims())

    # Create new dimensions in output file
    for dim in field_dims:
        if not any(x in dim.name for x in ('node', 'edge', 'face') +
                                          tuple(nc_file.dimensions.keys())):
            print('Adding previously undefined dimension {}'.format(dim.name))
            nc_file.createDimension(dim.name, dim.size)

    # Create new field variables in output file
    for field in fields:
        output_field_dims = []
        for dim in field.dimensions:
            if dim.endswith('face'):
                output_field_dims.append('nMesh2d_face')
            elif dim.endswith('edge'):
                output_field_dims.append('nMesh2d_edge')
            elif dim.endswith('node'):
                output_field_dims.append('nMesh2d_node')
            else:
                output_field_dims.append(dim)
        field_var = nc_file.createVariable(field.name, field.datatype, output_field_dims)
        field_var[:] = field[:]
        for att in field.ncattrs():
            if att == 'mesh':
                field_var.setncattr(att, output_mesh_name)
            elif att == 'coordinates':
                coord_names = field.getncattr(att).split(' ')
                new_coord_names = []
                for coord_name in coord_names:
                    if coord_name.startswith('Mesh2d'):
                        coord_name_comps = coord_name.split('_')
                        new_coord_names.append('_'.join([output_mesh_name] + coord_name_comps[-2:]))
                    else:
                        new_coord_names.append(coord_name)
                field_var.setncattr(att, ' '.join(new_coord_names))
            else:
                field_var.setncattr(att, field.getncattr(att))


def get_args():
    """
    Define and get command line arguments
    """

    parser = ap.ArgumentParser(description='Fix UGRID mesh description in LFRic output files')
    parser.add_argument('infile', help='LFRic UGRID input filename')
    parser.add_argument('meshfile', help='LFRic mesh filename')
    parser.add_argument('outfile', help='Output UGRID filename')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files')
    parser.add_argument('-q', '--quiet', action='store_true', help='Suppress output')
    return parser.parse_args()


def main():
    """
    Main routine
    """
    args = get_args()

    log_format = '[%(levelname)s] %(message)s'
    if args.quiet:
        logging.basicConfig(level=logging.WARNING, format=log_format)
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

    logging.info('LFRic UGRID input file: %s', args.infile)
    logging.info('LFRic mesh file: %s', args.meshfile)
    logging.info('Output file: %s', args.outfile)

    infile = nc.Dataset(args.infile, 'r')
    fields = get_fields(infile)
    vert_axes = get_vertical_axes(infile)
    time_vars = get_time_vars(infile)

    meshfile = nc.Dataset(args.meshfile, 'r')
    ugrid_vars = get_ugrid_mesh_description(meshfile)

    outfile = nc.Dataset(args.outfile, clobber=args.force, mode='w')
    write_file(outfile, ugrid_vars, vert_axes, time_vars, fields)

    outfile.close()
    meshfile.close()
    infile.close()


if __name__ == "__main__":
    main()
