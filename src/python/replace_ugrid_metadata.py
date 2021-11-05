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
    Scan input netCDF file and return a list of variables that are associated with model
    data fields, requiring certain netCDF variable attributes for identification
    """

    fields = []
    for var in nc_file.variables.values():
        if all(att in dir(var) for att in ['mesh', 'units', 'location']):
            fields.append(var)

    logging.info('Found %i field variables.', len(fields))

    return fields


def get_vertical_axes(nc_file):
    """
    Scan input netCDF file and return a list of vertical axis variables, requiring specific
    axis names
    """

    vertical_axes = []
    for var_name, var in nc_file.variables.items():
        if var_name in ('full_levels', 'half_levels'):
            vertical_axes.append(var)

    logging.info('Found %i vertical axes.', len(vertical_axes))

    return vertical_axes


def get_time_vars(nc_file):
    """
    Scan input netCDF file and return a list of time variables, requiring specific variable
    names
    """

    time_vars = []
    for var_name, var in nc_file.variables.items():
        if var_name in ('time_instant', 'time_instant_bounds'):
            time_vars.append(var)

    logging.info('Found %i time variables.', len(time_vars))

    return time_vars


def get_ugrid_mesh_description(nc_file):
    """
    Collect UGRID mesh description from mesh file and return a list of UGRID variables;
    only one mesh description is allowed
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

    # Extract topology dummy var from the first (and only) list element
    ugrid_topology_var = nc_file.variables[mesh_names[0]]

    # Collect expected UGRID location and connectivity variable names from the dummy var
    # attributes
    ugrid_varnames = set()
    ugrid_varnames.update(ugrid_topology_var.node_coordinates.split(" "))
    ugrid_varnames.update(ugrid_topology_var.face_coordinates.split(" "))
    ugrid_varnames.add(ugrid_topology_var.face_node_connectivity)
    ugrid_varnames.add(ugrid_topology_var.edge_node_connectivity)
    ugrid_varnames.add(ugrid_topology_var.face_edge_connectivity)
    ugrid_varnames.add(ugrid_topology_var.face_face_connectivity)

    # Assemble list of UGRID variables, including dummy var
    ugrid_vars = [ugrid_topology_var]
    for varname in ugrid_varnames:
        ugrid_vars.append(nc_file.variables[varname])

    return ugrid_vars


def write_ugrid_vars(nc_file, ugrid_vars, output_mesh_name):
    """
    Write UGRID dimension and variables
    """

    # Collect UGRID dimensions from variables
    ugrid_dims = set()
    for ugrid_var in ugrid_vars:
        ugrid_dims.update(ugrid_var.get_dims())

    mesh_name = ugrid_vars[0].name

    # Generate dictionary for translating between old and new UGRID dimension names
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


def generate_edge_coords(nc_file, output_mesh_name):
    """
    Generate edge coordinate field, if missing, by averaging over node locations
    """

    if 'edge_coordinates' not in dir(nc_file.variables[output_mesh_name]):
        logging.info('Generating edge coordinates...')
        edge_node_conn = nc_file.variables[output_mesh_name + '_edge_nodes']
        if 'start_index' in dir(edge_node_conn):
            start_idx = int(edge_node_conn.start_index)
        else:
            start_idx = 0
        node_lon = nc_file.variables[output_mesh_name + '_node_x']
        node_lat = nc_file.variables[output_mesh_name + '_node_y']
        edge_lon_var = nc_file.createVariable(output_mesh_name + '_edge_x', node_lon.datatype,
                                              'n' + output_mesh_name + '_edge')
        edge_lon_var[:] = 0.5*(node_lon[edge_node_conn[:,0] - start_idx] + \
                               node_lon[edge_node_conn[:,1] - start_idx])
        for att in node_lon.ncattrs():
            edge_lon_var.setncattr(att, node_lon.getncattr(att).replace('node', 'edge'))
        edge_lat_var = nc_file.createVariable(output_mesh_name + '_edge_y', node_lat.datatype,
                                              'n' + output_mesh_name + '_edge')
        edge_lat_var[:] = 0.5*(node_lat[edge_node_conn[:,0] - start_idx] + \
                               node_lat[edge_node_conn[:,1] - start_idx])
        for att in node_lat.ncattrs():
            edge_lat_var.setncattr(att, node_lat.getncattr(att).replace('node', 'edge'))
        nc_file.variables[output_mesh_name].setncattr('edge_coordinates', output_mesh_name + \
                                                      '_edge_x ' + output_mesh_name + '_edge_y')


def write_vars(nc_file, variables):
    """
    Write variable dimensions and data
    """

    dims = set()
    for var in variables:
        dims.update(var.get_dims())

    for dim in dims:
        nc_file.createDimension(dim.name, dim.size)

    for var in variables:
        out_var = nc_file.createVariable(var.name, var.datatype, var.dimensions)
        out_var[:] = var[:]
        for att in var.ncattrs():
            out_var.setncattr(att, var.getncattr(att))


def write_field_vars(nc_file, fields, output_mesh_name):
    """
    Write field dimensions and data
    """

    # Collect dimensions from input file
    field_dims = set()
    for field in fields:
        field_dims.update(field.get_dims())

    # Add any previously undefined dimensions - these should be field component dimensions,
    # but there is no way to be sure, so report them to the user
    for dim in field_dims:
        if not any(x in dim.name for x in ('node', 'edge', 'face') +
                                          tuple(nc_file.dimensions.keys())):
            msg = 'Adding previously undefined dimension {}'.format(dim.name)
            logging.info(msg)
            nc_file.createDimension(dim.name, dim.size)

    # Create new field variables in output file
    for field in fields:
        # Substitute dimension names
        output_field_dims = []
        for dim in field.dimensions:
            if dim.endswith('face'):
                output_field_dims.append('n' + output_mesh_name + '_face')
            elif dim.endswith('edge'):
                output_field_dims.append('n' + output_mesh_name + '_edge')
            elif dim.endswith('node'):
                output_field_dims.append('n' + output_mesh_name + '_node')
            else:
                output_field_dims.append(dim)

        field_var = nc_file.createVariable(field.name, field.datatype, output_field_dims)
        field_var[:] = field[:]

        for att in field.ncattrs():
            # Substitute old UGRID mesh name
            if att == 'mesh':
                field_var.setncattr(att, output_mesh_name)
            # Reconstruct coordinate names with new UGRID mesh name, if applicable
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
    parser.add_argument('-m', '--meshname', default='Mesh2d', help='Output UGRID mesh name')
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

    # Load LFRic outut file and retrieve data fields, vertical axes, and time variables
    infile = nc.Dataset(args.infile, 'r')
    fields = get_fields(infile)
    vert_axes = get_vertical_axes(infile)
    time_vars = get_time_vars(infile)

    # Load LFRic mesh file and retrieve UGRID variables
    meshfile = nc.Dataset(args.meshfile, 'r')
    ugrid_vars = get_ugrid_mesh_description(meshfile)

    # Open output file and write variables and dimensions
    outfile = nc.Dataset(args.outfile, clobber=args.force, mode='w')
    write_ugrid_vars(outfile, ugrid_vars, args.meshname)
    generate_edge_coords(outfile, args.meshname)
    write_vars(outfile, vert_axes)
    write_vars(outfile, time_vars)
    write_field_vars(outfile, fields, args.meshname)

    outfile.close()
    meshfile.close()
    infile.close()


if __name__ == "__main__":
    main()
