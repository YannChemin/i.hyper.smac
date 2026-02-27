"""
Name:       i.hyper.smac test
Purpose:    Tests i.hyper.smac atmospheric correction module.

Author:     Yann Chemin - yann.chemin@gmail.com
Copyright:  (C) 2025 by Yann Chemin and the GRASS Development Team
Licence:    This program is free software under the GNU General Public
            License (>=v2). Read the file COPYING that comes with GRASS
            for details.
"""

from grass.gunittest.case import TestCase
from grass.gunittest.gutils import is_map_in_mapset
from grass.script import array as garray
import grass.script as gs
import numpy as np
import os
import tempfile


class TestHyperSmacBasic(TestCase):
    """Test basic functionality of i.hyper.smac"""

    @classmethod
    def setUpClass(cls):
        """Set up test data and temporary region"""
        # Create temporary region
        cls.runModule("g.region", s=0, n=100, w=0, e=100, res=1, flags="p")
        cls.use_temp_region()
        
        # Create test DEM
        cls.runModule("r.mapcalc", expression="test_dem = 100.0 + row() * 0.1")
        
        # Create synthetic hyperspectral data for testing
        # Simulate a simple 3-band hyperspectral cube
        cls.runModule("r.mapcalc", expression="test_band_400 = 0.1 * sin(row() * 0.1) + 0.2")
        cls.runModule("r.mapcalc", expression="test_band_550 = 0.15 * sin(row() * 0.1) + 0.25")
        cls.runModule("r.mapcalc", expression="test_band_700 = 0.12 * sin(row() * 0.1) + 0.18")
        
        # Create 3D raster from individual bands
        cls.runModule("r.in.xyz", 
                     input="-", 
                     output="test_hyper",
                     separator="comma",
                     x="1", y="1", z="1",
                     flags="t",
                     stdin_="1,1,0.15\n2,2,0.20\n3,3,0.18\n")

    @classmethod
    def tearDownClass(cls):
        """Clean up test data"""
        maps_to_remove = [
            "test_dem", "test_band_400", "test_band_550", "test_band_700",
            "test_hyper", "test_corrected"
        ]
        cls.runModule(
            "g.remove",
            flags="f",
            type="raster",
            name=",".join(maps_to_remove),
        )
        cls.del_temp_region()

    def test_module_exists(self):
        """Test that i.hyper.smac module exists and can be called"""
        self.assertModule("i.hyper.smac", flags="h")

    def test_required_parameters(self):
        """Test that required parameters are properly validated"""
        # Test with missing input
        self.assertModuleFail("i.hyper.smac", dem="test_dem", solar_zenith=45)
        
        # Test with missing DEM
        self.assertModuleFail("i.hyper.smac", input="test_hyper", solar_zenith=45)
        
        # Test with missing solar zenith
        self.assertModuleFail("i.hyper.smac", input="test_hyper", dem="test_dem")

    def test_basic_atmospheric_correction(self):
        """Test basic atmospheric correction functionality"""
        # Run atmospheric correction with minimal parameters
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper",
            output="test_corrected",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            water_vapor=1.5,
            ozone=0.3
        )
        
        # Check that output was created
        self.assertTrue(is_map_in_mapset("test_corrected", type="raster3d"))

    def test_parameter_ranges(self):
        """Test parameter validation and ranges"""
        # Test invalid solar zenith
        self.assertModuleFail(
            "i.hyper.smac",
            input="test_hyper",
            output="test_corrected",
            dem="test_dem",
            solar_zenith=95,  # Invalid: > 90 degrees
            aod=0.2
        )
        
        # Test invalid AOD
        self.assertModuleFail(
            "i.hyper.smac",
            input="test_hyper", 
            output="test_corrected",
            dem="test_dem",
            solar_zenith=45,
            aod=-0.1  # Invalid: negative AOD
        )

    def test_smart_lut_parameter(self):
        """Test smart LUT functionality"""
        # Test with smart LUT enabled
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper",
            output="test_corrected",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            smart_lut="yes"
        )
        
        # Test with smart LUT disabled
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper",
            output="test_corrected2",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            smart_lut="no"
        )

    def test_parallel_lut_parameter(self):
        """Test parallel LUT functionality"""
        # Test with parallel LUT enabled
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper",
            output="test_corrected",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            parallel_lut="enabled"
        )
        
        # Test with parallel LUT disabled
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper",
            output="test_corrected2",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            parallel_lut="disabled"
        )


class TestHyperSmacAODEstimation(TestCase):
    """Test AOD estimation functionality"""

    @classmethod
    def setUpClass(cls):
        """Set up test data for AOD estimation"""
        cls.runModule("g.region", s=0, n=50, w=0, e=50, res=1, flags="p")
        cls.use_temp_region()
        
        # Create test DEM
        cls.runModule("r.mapcalc", expression="test_dem = 120.0")
        
        # Create synthetic hyperspectral data suitable for AOD estimation
        # Simulate dark target pixels
        cls.runModule("r.mapcalc", expression="blue_band = 0.05 + rand(0, 0.02)")
        cls.runModule("r.mapcalc", expression="red_band = 0.08 + rand(0, 0.02)")
        cls.runModule("r.mapcalc", expression="nir_band = 0.15 + rand(0, 0.02)")
        cls.runModule("r.mapcalc", expression="swir_band = 0.12 + rand(0, 0.02)")
        
        # Create 3D raster
        cls.runModule("r.in.xyz",
                     input="-",
                     output="test_hyper_aod", 
                     separator="comma",
                     x="1", y="1", z="1",
                     flags="t",
                     stdin_="1,1,0.08\n2,2,0.10\n3,3,0.12\n")

    @classmethod
    def tearDownClass(cls):
        """Clean up test data"""
        maps_to_remove = [
            "test_dem", "blue_band", "red_band", "nir_band", "swir_band",
            "test_hyper_aod", "test_corrected_aod"
        ]
        cls.runModule(
            "g.remove",
            flags="f",
            type="raster",
            name=",".join(maps_to_remove),
        )
        cls.del_temp_region()

    def test_aod_estimation(self):
        """Test automatic AOD estimation"""
        # Run without providing AOD - should estimate automatically
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_aod",
            output="test_corrected_aod",
            dem="test_dem",
            solar_zenith=45,
            water_vapor=1.5,
            ozone=0.3
        )
        
        # Check that output was created
        self.assertTrue(is_map_in_mapset("test_corrected_aod", type="raster3d"))

    def test_aod_map_output(self):
        """Test AOD map generation"""
        # Run with AOD estimation
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_aod",
            output="test_corrected_aod",
            dem="test_dem",
            solar_zenith=45,
            aod_map="test_aod_map"
        )
        
        # Check that AOD map was created
        self.assertTrue(is_map_in_mapset("test_aod_map", type="raster"))


class TestHyperSmacWaterVapor(TestCase):
    """Test water vapor estimation functionality"""

    @classmethod
    def setUpClass(cls):
        """Set up test data for water vapor estimation"""
        cls.runModule("g.region", s=0, n=30, w=0, e=30, res=1, flags="p")
        cls.use_temp_region()
        
        # Create test DEM
        cls.runModule("r.mapcalc", expression="test_dem = 95.0")
        
        # Create synthetic data with water vapor absorption features
        cls.runModule("r.mapcalc", expression="band_940 = 0.10 * exp(-row() * 0.05)")
        cls.runModule("r.mapcalc", expression="band_1130 = 0.12 * exp(-row() * 0.03)")
        
        # Create 3D raster
        cls.runModule("r.in.xyz",
                     input="-",
                     output="test_hyper_wv",
                     separator="comma", 
                     x="1", y="1", z="1",
                     flags="t",
                     stdin_="1,1,0.11\n2,2,0.13\n3,3,0.09\n")

    @classmethod
    def tearDownClass(cls):
        """Clean up test data"""
        maps_to_remove = [
            "test_dem", "band_940", "band_1130", "test_hyper_wv",
            "test_corrected_wv"
        ]
        cls.runModule(
            "g.remove",
            flags="f",
            type="raster",
            name=",".join(maps_to_remove),
        )
        cls.del_temp_region()

    def test_water_vapor_estimation(self):
        """Test water vapor estimation methods"""
        # Test joint estimation method
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_wv",
            output="test_corrected_wv",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            water_vapor="joint"
        )
        
        # Test 940nm method
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_wv",
            output="test_corrected_wv2",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            water_vapor="940nm"
        )

    def test_wvc_map_output(self):
        """Test water vapor content map generation"""
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_wv",
            output="test_corrected_wv",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            wvc_map="test_wvc_map"
        )
        
        # Check that WVC map was created
        self.assertTrue(is_map_in_mapset("test_wvc_map", type="raster"))


class TestHyperSmacAdvanced(TestCase):
    """Test advanced functionality and edge cases"""

    @classmethod
    def setUpClass(cls):
        """Set up test data for advanced tests"""
        cls.runModule("g.region", s=0, n=20, w=0, e=20, res=1, flags="p")
        cls.use_temp_region()
        
        # Create test DEM with elevation variation
        cls.runModule("r.mapcalc", expression="test_dem = 100.0 + col() * 0.5")
        
        # Create simple test hyperspectral data
        cls.runModule("r.in.xyz",
                     input="-",
                     output="test_hyper_adv",
                     separator="comma",
                     x="1", y="1", z="1", 
                     flags="t",
                     stdin_="1,1,0.15\n2,2,0.18\n3,3,0.12\n")

    @classmethod
    def tearDownClass(cls):
        """Clean up test data"""
        maps_to_remove = [
            "test_dem", "test_hyper_adv", "test_corrected_adv",
            "test_uncertainty"
        ]
        cls.runModule(
            "g.remove",
            flags="f",
            type="raster",
            name=",".join(maps_to_remove),
        )
        cls.del_temp_region()

    def test_uncertainty_estimation(self):
        """Test uncertainty estimation functionality"""
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_adv",
            output="test_corrected_adv",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            compute_uncertainty=True,
            output_uncertainty="test_uncertainty"
        )
        
        # Check that uncertainty map was created
        self.assertTrue(is_map_in_mapset("test_uncertainty", type="raster3d"))

    def test_adjacency_correction(self):
        """Test adjacency correction functionality"""
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_adv",
            output="test_corrected_adv",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            adjacency_psf=1.0
        )

    def test_aerosol_models(self):
        """Test different aerosol models"""
        aerosol_models = ["continental", "maritime", "urban", "desert"]
        
        for model in aerosol_models:
            self.assertModule(
                "i.hyper.smac",
                input="test_hyper_adv",
                output=f"test_corrected_{model}",
                dem="test_dem",
                solar_zenith=45,
                aod=0.2,
                aerosol_model=model
            )

    def test_opencl_functionality(self):
        """Test OpenCL acceleration"""
        # Test with OpenCL enabled
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_adv",
            output="test_corrected_opencl",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            opencl_device="cpu",
            opencl_memory=512
        )

    def test_force_regenerate(self):
        """Test cache regeneration functionality"""
        # First run
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_adv",
            output="test_corrected_adv",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2
        )
        
        # Second run with force regeneration
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_adv",
            output="test_corrected_adv2",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            clear_lut_cache=True
        )


class TestHyperSmacPerformance(TestCase):
    """Test performance and resource usage"""

    @classmethod
    def setUpClass(cls):
        """Set up test data for performance tests"""
        cls.runModule("g.region", s=0, n=10, w=0, e=10, res=1, flags="p")
        cls.use_temp_region()
        
        # Create small test dataset
        cls.runModule("r.mapcalc", expression="test_dem = 100.0")
        cls.runModule("r.in.xyz",
                     input="-",
                     output="test_hyper_perf",
                     separator="comma",
                     x="1", y="1", z="1",
                     flags="t",
                     stdin_="1,1,0.15\n2,2,0.18\n")

    @classmethod
    def tearDownClass(cls):
        """Clean up test data"""
        maps_to_remove = ["test_dem", "test_hyper_perf", "test_corrected_perf"]
        cls.runModule(
            "g.remove",
            flags="f",
            type="raster",
            name=",".join(maps_to_remove),
        )
        cls.del_temp_region()

    def test_smart_vs_standard_lut(self):
        """Test performance difference between smart and standard LUT"""
        import time
        
        # Test standard LUT
        start_time = time.time()
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_perf",
            output="test_corrected_standard",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            smart_lut="no"
        )
        standard_time = time.time() - start_time
        
        # Test smart LUT
        start_time = time.time()
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_perf",
            output="test_corrected_smart",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            smart_lut="yes"
        )
        smart_time = time.time() - start_time
        
        # Smart LUT should generally be faster (though this is a small test)
        self.assertLess(smart_time, standard_time * 2)  # Allow some variance

    def test_parallel_vs_serial_lut(self):
        """Test performance difference between parallel and serial LUT"""
        import time
        
        # Test serial LUT
        start_time = time.time()
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_perf",
            output="test_corrected_serial",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            parallel_lut="disabled"
        )
        serial_time = time.time() - start_time
        
        # Test parallel LUT
        start_time = time.time()
        self.assertModule(
            "i.hyper.smac",
            input="test_hyper_perf",
            output="test_corrected_parallel",
            dem="test_dem",
            solar_zenith=45,
            aod=0.2,
            parallel_lut="enabled"
        )
        parallel_time = time.time() - start_time
        
        # Parallel should generally be faster on multi-core systems
        # (though small test may not show significant difference)
        self.assertLess(parallel_time, serial_time * 2)


if __name__ == "__main__":
    from grass.gunittest.main import test
    test()
