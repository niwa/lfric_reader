"""
Unit tests for script replace_ugrid_metadata
"""
import os
import sys
import unittest
import netCDF4 as nc

# Make parent directory available for import
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

import replace_ugrid_metadata as rum

class ReplaceUGRIDMetadataUnitTest(unittest.TestCase):
    """
    Unit tests for functions in replace_ugrid_metadata
    """

    def test_get_fields(self):
        """
        Unit test for get_fields function
        """

        nc_file = nc.Dataset('test_get_fields_nc', clobber=False, mode='w', diskless=True,
                             persist=False)

        # Valid variable
        nc_var = nc_file.createVariable('var1', 'i')
        nc_var.setncattr('mesh', '123')
        nc_var.setncattr('units', '123')
        nc_var.setncattr('location', '123')

        # Invalid variable
        nc_file.createVariable('var2', 'i')

        fields = rum.get_fields(nc_file)

        self.assertIsInstance(fields, list)
        self.assertEqual(len(fields), 1)
        self.assertEqual(fields[0].name, 'var1')

        nc_file.close()


    def test_get_vertical_axes(self):
        """
        Unit test for get_vertical_axes function
        """

        nc_file = nc.Dataset('test_get_vertical_axes.nc', clobber=False, mode='w', diskless=True,
                             persist=False)

        # Valid variables
        nc_file.createVariable('half_levels', 'i')
        nc_file.createVariable('full_levels', 'i')

        # Invalid variable
        nc_file.createVariable('var3', 'i')

        axes = rum.get_vertical_axes(nc_file)

        self.assertIsInstance(axes, list)
        self.assertEqual(len(axes), 2)
        self.assertEqual(axes[0].name, 'half_levels')
        self.assertEqual(axes[1].name, 'full_levels')

        nc_file.close()


    def test_get_time_vars(self):
        """
        Unit test for get_time_vars function
        """

        nc_file = nc.Dataset('test_get_time_vars.nc', clobber=False, mode='w', diskless=True,
                             persist=False)

        # Valid variables
        nc_file.createVariable('time_instant', 'i')
        nc_file.createVariable('time_instant_bounds', 'i')

        # Invalid variable
        nc_file.createVariable('var3', 'i')

        axes = rum.get_time_vars(nc_file)

        self.assertIsInstance(axes, list)
        self.assertEqual(len(axes), 2)
        self.assertEqual(axes[0].name, 'time_instant')
        self.assertEqual(axes[1].name, 'time_instant_bounds')

        nc_file.close()


    def test_get_ugrid_mesh_description(self):
        """
        Unit test for get_ugrid_mesh_description
        """

        nc_file = nc.Dataset('test_get_ugrid_1.nc', clobber=False, mode='w', diskless=True,
                             persist=False)

        # Only one mesh topology allowed per file
        nc_var1 = nc_file.createVariable('var1', 'i')
        nc_var1.setncattr('cf_role', 'mesh_topology')
        nc_var2 = nc_file.createVariable('var2', 'i')
        nc_var2.setncattr('cf_role', 'mesh_topology')
        self.assertRaises(RuntimeError, rum.get_ugrid_mesh_description, nc_file)

        nc_file.close()

        nc_file = nc.Dataset('test_get_ugrid_2.nc', clobber=False, mode='w', diskless=True,
                             persist=False)

        # Set up valid UGRID variables
        nc_var1 = nc_file.createVariable('var1', 'i')
        att_dict = {
            'cf_role' : 'mesh_topology',
            'node_coordinates' : 'node_x node_y',
            'face_coordinates' : 'face_x face_y',
            'face_node_connectivity' : 'face_nodes',
            'edge_node_connectivity' : 'edge_nodes',
            'face_edge_connectivity' : 'face_edges',
            'face_face_connectivity' : 'face_links',
        }
        for key, value in att_dict.items():
            nc_var1.setncattr(key, value)
        var_list = ['node_x', 'node_y', 'face_x', 'face_y', 'face_nodes',
                    'edge_nodes', 'face_edges', 'face_links', 'extra_var']
        for var in var_list:
            nc_file.createVariable(var, 'i')

        ugrid = rum.get_ugrid_mesh_description(nc_file)

        self.assertIsInstance(ugrid, list)
        self.assertEqual(len(ugrid), 9)
        for nc_var in ugrid:
            self.assertIn(nc_var.name, ['var1'] + var_list[:-1])

        nc_file.close()


    def test_write_ugrid_vars(self):
        """
        Unit test for write_ugrid_vars
        """

        # Define a few source dims variables
        nc_file_src = nc.Dataset('test_write_ugrid_src.nc', clobber=False, mode='w', diskless=True,
                                 persist=False)
        nc_file_src.createDimension('ugrid_mesh_dim1', 1)
        nc_file_src.createDimension('ugrid_mesh_dim2', 2)
        nc_file_src.createDimension('ugrid_mesh_dim3', 3)
        var_list = []
        var_list.append(nc_file_src.createVariable('ugrid_mesh', 'i'))
        var_list.append(nc_file_src.createVariable('var1', 'i', 'ugrid_mesh_dim1'))
        var_list.append(nc_file_src.createVariable('var2', 'i', 'ugrid_mesh_dim2'))
        var_list[-1][:] = 2
        var_list[-1].setncattr('abc', 'ugrid_mesh_att')
        var_list[-1].setncattr('def', 'ijk')

        nc_file_dst = nc.Dataset('test_write_ugrid_dst.nc', clobber=False, mode='w', diskless=True,
                                 persist=False)

        rum.write_ugrid_vars(nc_file_dst, var_list, 'new_mesh')

        self.assertEqual(len(nc_file_dst.dimensions), 2)
        self.assertIn('new_mesh_dim1', nc_file_dst.dimensions)
        self.assertIn('new_mesh_dim2', nc_file_dst.dimensions)
        self.assertNotIn('new_mesh_dim3', nc_file_dst.dimensions)
        self.assertEqual(len(nc_file_dst.variables), 3)
        self.assertIn('new_mesh', nc_file_dst.variables)
        self.assertIn('var1', nc_file_dst.variables)
        self.assertIn('var2', nc_file_dst.variables)
        self.assertEqual(nc_file_dst.variables['var2'][:].size,
                         nc_file_dst.dimensions['new_mesh_dim2'].size)
        self.assertListEqual(list(nc_file_dst.variables['var2'][:]), [2, 2])
        self.assertEqual(nc_file_dst.variables['var2'].getncattr('abc'), 'new_mesh_att')
        self.assertEqual(nc_file_dst.variables['var2'].getncattr('def'), 'ijk')

        nc_file_src.close()
        nc_file_dst.close()


    def test_generate_edge_coords(self):
        """
        Unit test for generate_edge_coords
        """

        nc_file = nc.Dataset('test_edge_coords_1.nc', clobber=False, mode='w', diskless=True,
                             persist=False)
        nc_var = nc_file.createVariable('ugrid_mesh', 'i')
        nc_var.setncattr('edge_coordinates', 'var1 var2')
        nc_file.createVariable('var1', 'i')
        nc_file.createVariable('var2', 'i')

        rum.generate_edge_coords(nc_file, 'ugrid_mesh')

        # File is not modified if UGRID mesh has edge_coordinates attribute
        self.assertEqual(len(nc_file.variables), 3)
        self.assertIn('ugrid_mesh', nc_file.variables)
        self.assertEqual(nc_file.variables['ugrid_mesh'].getncattr('edge_coordinates'), 'var1 var2')
        self.assertIn('var1', nc_file.variables)
        self.assertIn('var2', nc_file.variables)

        nc_file.close()

        # Generate file with a single edge
        nc_file = nc.Dataset('test_edge_coords_2.nc', clobber=False, mode='w', diskless=True,
                             persist=False)
        nc_file.createDimension('Two', 2)
        nc_file.createDimension('nugrid_mesh_node', 2)
        nc_file.createDimension('nugrid_mesh_edge', 1)
        nc_var = nc_file.createVariable('ugrid_mesh', 'i')
        nc_var.setncattr('node_coordinates', 'ugrid_mesh_node_x ugrid_mesh_node_y')
        nc_var.setncattr('edge_node_connectivity', 'ugrid_mesh_edge_nodes')
        node_x_var = nc_file.createVariable('ugrid_mesh_node_x', 'f', 'nugrid_mesh_node')
        node_x_var[:] = [0.0, 1.0]
        node_y_var = nc_file.createVariable('ugrid_mesh_node_y', 'f', 'nugrid_mesh_node')
        node_y_var[:] = [0.0, 0.0]
        en_conn_var = nc_file.createVariable('ugrid_mesh_edge_nodes', 'i',
                                             ('nugrid_mesh_edge', 'Two'))
        en_conn_var[:] = [0, 1]

        rum.generate_edge_coords(nc_file, 'ugrid_mesh')

        # File has additional edge coordinate variable
        self.assertEqual(len(nc_file.variables), 6)
        self.assertEqual(nc_file.variables['ugrid_mesh'].getncattr('edge_coordinates'),
                         'ugrid_mesh_edge_x ugrid_mesh_edge_y')
        self.assertIn('ugrid_mesh_edge_x', nc_file.variables)
        self.assertIn('ugrid_mesh_edge_y', nc_file.variables)
        self.assertAlmostEqual(nc_file.variables['ugrid_mesh_edge_x'][:], 0.5)
        self.assertAlmostEqual(nc_file.variables['ugrid_mesh_edge_y'][:], 0.0)


    def test_write_vars(self):
        """
        Unit test for write_vars
        """

        # Define a few source dims and variables
        nc_file_src = nc.Dataset('test_write_vars_src.nc', clobber=False, mode='w', diskless=True,
                                 persist=False)
        nc_file_src.createDimension('dim1', 1)
        nc_file_src.createDimension('dim2', 2)
        nc_file_src.createDimension('dim3', 3)
        var_list = []
        var_list.append(nc_file_src.createVariable('var1', 'i', 'dim1'))
        var_list.append(nc_file_src.createVariable('var2', 'i', 'dim2'))
        var_list.append(nc_file_src.createVariable('var3', 'i'))
        var_list[0][:] = 1
        var_list[-1].setncattr('abc', 'def')


        nc_file_dst = nc.Dataset('test_write_vars_dst.nc', clobber=False, mode='w', diskless=True,
                                 persist=False)

        rum.write_vars(nc_file_dst, var_list)

        self.assertEqual(len(nc_file_dst.dimensions), 2)
        self.assertIn('dim1', nc_file_dst.dimensions)
        self.assertIn('dim2', nc_file_dst.dimensions)
        self.assertNotIn('dim3', nc_file_dst.dimensions)
        self.assertEqual(len(nc_file_dst.variables), 3)
        self.assertIn('var1', nc_file_dst.variables)
        self.assertIn('var2', nc_file_dst.variables)
        self.assertIn('var3', nc_file_dst.variables)
        self.assertEqual(nc_file_dst.variables['var1'][:].size, nc_file_dst.dimensions['dim1'].size)
        self.assertTrue(nc_file_dst.variables['var1'][:] == 1)
        self.assertEqual(nc_file_dst.variables['var3'].getncattr('abc'), 'def')

        nc_file_src.close()
        nc_file_dst.close()


    def test_write_field_vars(self):
        """
        Unit test for write_field_vars
        """

        nc_file_src = nc.Dataset('test_write_field_vars_src.nc', clobber=False, mode='w',
                                 diskless=True, persist=False)
        nc_file_src.createDimension('node', 4)
        nc_file_src.createDimension('edge', 4)
        nc_file_src.createDimension('face', 1)
        nc_file_src.createDimension('component_dim', 1)
        var_list = []
        var_list.append(nc_file_src.createVariable('node_var', 'f', 'node'))
        var_list.append(nc_file_src.createVariable('edge', 'f', 'edge'))
        var_list.append(nc_file_src.createVariable('face_var', 'f', ('face', 'component_dim')))
        var_list[-1].setncattr('mesh', 'Mesh2d')
        var_list[-1].setncattr('coordinates', 'Mesh2d_face_x Mesh2d_face_y')

        nc_file_dst = nc.Dataset('test_write_field_vars_dst.nc', clobber=False, mode='w',
                                 diskless=True, persist=False)
        nc_file_dst.createDimension('nugrid_mesh_node', 4)
        nc_file_dst.createDimension('nugrid_mesh_edge', 4)
        nc_file_dst.createDimension('nugrid_mesh_face', 1)

        rum.write_field_vars(nc_file_dst, var_list, 'ugrid_mesh')

        self.assertEqual(len(nc_file_dst.dimensions), 4)
        self.assertIn('component_dim', nc_file_dst.dimensions)
        self.assertEqual(nc_file_dst['face_var'].getncattr('mesh'), 'ugrid_mesh')
        self.assertEqual(nc_file_dst['face_var'].getncattr('coordinates'),
                         'ugrid_mesh_face_x ugrid_mesh_face_y')

        nc_file_src.close()
        nc_file_dst.close()
