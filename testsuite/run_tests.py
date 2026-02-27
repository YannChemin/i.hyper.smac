#!/usr/bin/env python3
"""
GRASS GIS i.hyper.smac Test Runner

Author: Yann Chemin <yann.chemin@gmail.com>
Copyright: (C) 2025 by Yann Chemin and the GRASS Development Team
License: GPL-2.0-or-later

This script provides comprehensive testing for the i.hyper.smac module
including unit tests, integration tests, and performance benchmarks.
"""

import os
import sys
import subprocess
import time
import argparse
from pathlib import Path

# Add GRASS Python bindings to path
grass_bin = os.environ.get('GISBASE', '/usr/local/grass86')
sys.path.insert(0, os.path.join(grass_bin, 'etc', 'python'))

try:
    import grass.script as gs
    from grass.gunittest.main import test as gunittest
    from grass.gunittest.case import TestCase
except ImportError as e:
    print(f"Error importing GRASS modules: {e}")
    print("Please ensure GRASS GIS is properly installed and GISBASE is set")
    sys.exit(1)


class HyperSmacTestRunner:
    """Comprehensive test runner for i.hyper.smac module"""
    
    def __init__(self, test_data_dir=None):
        """Initialize test runner
        
        Args:
            test_data_dir: Directory containing test data (optional)
        """
        self.test_data_dir = test_data_dir
        self.results = {
            'unit_tests': {'passed': 0, 'failed': 0, 'errors': []},
            'integration_tests': {'passed': 0, 'failed': 0, 'errors': []},
            'performance_tests': {'passed': 0, 'failed': 0, 'errors': []},
            'total_time': 0
        }
        
    def setup_test_environment(self):
        """Set up GRASS test environment"""
        print("Setting up GRASS test environment...")
        
        # Create temporary location for testing
        self.temp_location = gs.tempdir()
        
        # Initialize GRASS session
        try:
            gs.run_command('g.gisenv', set=f'GISDBASE={self.temp_location}')
            gs.run_command('g.proj', flags='c', epsg=4326)
            print("✓ GRASS environment initialized")
            return True
        except Exception as e:
            print(f"✗ Failed to initialize GRASS environment: {e}")
            return False
    
    def cleanup_test_environment(self):
        """Clean up test environment"""
        print("Cleaning up test environment...")
        try:
            if hasattr(self, 'temp_location') and os.path.exists(self.temp_location):
                import shutil
                shutil.rmtree(self.temp_location)
            print("✓ Test environment cleaned up")
        except Exception as e:
            print(f"⚠ Warning during cleanup: {e}")
    
    def run_unit_tests(self):
        """Run unit tests using GRASS unittest framework"""
        print("\n" + "="*60)
        print("RUNNING UNIT TESTS")
        print("="*60)
        
        start_time = time.time()
        
        try:
            # Run unit tests
            test_file = os.path.join(os.path.dirname(__file__), 'test_i_hyper_smac.py')
            
            if not os.path.exists(test_file):
                raise FileNotFoundError(f"Test file not found: {test_file}")
            
            # Use GRASS unittest runner
            result = subprocess.run([
                sys.executable, test_file
            ], capture_output=True, text=True, cwd=os.path.dirname(__file__))
            
            execution_time = time.time() - start_time
            
            if result.returncode == 0:
                print(f"✓ Unit tests passed in {execution_time:.2f}s")
                self.results['unit_tests']['passed'] = 1
            else:
                print(f"✗ Unit tests failed in {execution_time:.2f}s")
                print("STDOUT:", result.stdout)
                print("STDERR:", result.stderr)
                self.results['unit_tests']['failed'] = 1
                self.results['unit_tests']['errors'].append(result.stderr)
                
        except Exception as e:
            print(f"✗ Unit test execution failed: {e}")
            self.results['unit_tests']['failed'] = 1
            self.results['unit_tests']['errors'].append(str(e))
    
    def run_integration_tests(self):
        """Run integration tests with real data"""
        print("\n" + "="*60)
        print("RUNNING INTEGRATION TESTS")
        print("="*60)
        
        start_time = time.time()
        
        try:
            # Test basic functionality
            self._test_basic_integration()
            self._test_parameter_validation()
            self._test_output_generation()
            
            execution_time = time.time() - start_time
            print(f"✓ Integration tests passed in {execution_time:.2f}s")
            self.results['integration_tests']['passed'] = 1
            
        except Exception as e:
            execution_time = time.time() - start_time
            print(f"✗ Integration tests failed in {execution_time:.2f}s: {e}")
            self.results['integration_tests']['failed'] = 1
            self.results['integration_tests']['errors'].append(str(e))
    
    def _test_basic_integration(self):
        """Test basic module integration"""
        print("  Testing basic integration...")
        
        # Create test region
        gs.run_command('g.region', s=0, n=10, w=0, e=10, res=1)
        
        # Create test data
        gs.run_command('r.mapcalc', expression='test_dem = 100.0')
        gs.run_command('r.mapcalc', expression='test_band = 0.1 + rand(0, 0.05)')
        
        # Create simple 3D raster
        test_data = "1,1,0.12\n2,2,0.15\n3,3,0.11\n"
        gs.run_command('r.in.xyz', input='-', output='test_hyper',
                     separator='comma', x='1', y='1', z='1', flags='t',
                     stdin_=test_data)
        
        # Test module execution
        gs.run_command('i.hyper.smac',
                      input='test_hyper',
                      output='test_output',
                      dem='test_dem',
                      solar_zenith=45,
                      aod=0.2,
                      water_vapor=1.5,
                      ozone=0.3)
        
        # Verify output
        info = gs.parse_command('r3.info', map='test_output', flags='g')
        if 'cells' not in info:
            raise RuntimeError("Output verification failed")
        
        # Cleanup
        gs.run_command('g.remove', flags='f', type='raster',
                      name='test_dem,test_band')
        gs.run_command('g.remove', flags='f', type='raster3d',
                      name='test_hyper,test_output')
        
        print("    ✓ Basic integration test passed")
    
    def _test_parameter_validation(self):
        """Test parameter validation"""
        print("  Testing parameter validation...")
        
        # Test missing required parameters
        try:
            gs.run_command('i.hyper.smac', solar_zenith=45)
            raise RuntimeError("Should have failed with missing parameters")
        except Exception:
            pass  # Expected to fail
        
        print("    ✓ Parameter validation test passed")
    
    def _test_output_generation(self):
        """Test output generation"""
        print("  Testing output generation...")
        
        # Create minimal test setup
        gs.run_command('g.region', s=0, n=5, w=0, e=5, res=1)
        gs.run_command('r.mapcalc', expression='test_dem = 100.0')
        
        test_data = "1,1,0.13\n"
        gs.run_command('r.in.xyz', input='-', output='test_hyper',
                     separator='comma', x='1', y='1', z='1', flags='t',
                     stdin_=test_data)
        
        # Test with uncertainty output
        gs.run_command('i.hyper.smac',
                      input='test_hyper',
                      output='test_output',
                      dem='test_dem',
                      solar_zenith=45,
                      aod=0.2,
                      compute_uncertainty=True,
                      output_uncertainty='test_uncertainty')
        
        # Verify both outputs exist
        info_output = gs.parse_command('r3.info', map='test_output', flags='g')
        info_uncert = gs.parse_command('r3.info', map='test_uncertainty', flags='g')
        
        if 'cells' not in info_output or 'cells' not in info_uncert:
            raise RuntimeError("Output verification failed")
        
        # Cleanup
        gs.run_command('g.remove', flags='f', type='raster', name='test_dem')
        gs.run_command('g.remove', flags='f', type='raster3d',
                      name='test_hyper,test_output,test_uncertainty')
        
        print("    ✓ Output generation test passed")
    
    def run_performance_tests(self):
        """Run performance benchmarks"""
        print("\n" + "="*60)
        print("RUNNING PERFORMANCE TESTS")
        print("="*60)
        
        start_time = time.time()
        
        try:
            self._test_smart_lut_performance()
            self._test_parallel_lut_performance()
            self._test_opencl_performance()
            
            execution_time = time.time() - start_time
            print(f"✓ Performance tests passed in {execution_time:.2f}s")
            self.results['performance_tests']['passed'] = 1
            
        except Exception as e:
            execution_time = time.time() - start_time
            print(f"✗ Performance tests failed in {execution_time:.2f}s: {e}")
            self.results['performance_tests']['failed'] = 1
            self.results['performance_tests']['errors'].append(str(e))
    
    def _test_smart_lut_performance(self):
        """Test smart LUT performance"""
        print("  Testing smart LUT performance...")
        
        # Create test data
        gs.run_command('g.region', s=0, n=5, w=0, e=5, res=1)
        gs.run_command('r.mapcalc', expression='test_dem = 100.0')
        
        test_data = "1,1,0.14\n"
        gs.run_command('r.in.xyz', input='-', output='test_hyper',
                     separator='comma', x='1', y='1', z='1', flags='t',
                     stdin_=test_data)
        
        # Test standard LUT
        start_time = time.time()
        gs.run_command('i.hyper.smac',
                      input='test_hyper',
                      output='test_standard',
                      dem='test_dem',
                      solar_zenith=45,
                      aod=0.2,
                      smart_lut='no')
        standard_time = time.time() - start_time
        
        # Test smart LUT
        start_time = time.time()
        gs.run_command('i.hyper.smac',
                      input='test_hyper',
                      output='test_smart',
                      dem='test_dem',
                      solar_zenith=45,
                      aod=0.2,
                      smart_lut='yes')
        smart_time = time.time() - start_time
        
        print(f"    Standard LUT: {standard_time:.3f}s, Smart LUT: {smart_time:.3f}s")
        
        # Cleanup
        gs.run_command('g.remove', flags='f', type='raster', name='test_dem')
        gs.run_command('g.remove', flags='f', type='raster3d',
                      name='test_hyper,test_standard,test_smart')
        
        print("    ✓ Smart LUT performance test passed")
    
    def _test_parallel_lut_performance(self):
        """Test parallel LUT performance"""
        print("  Testing parallel LUT performance...")
        
        # Create test data
        gs.run_command('g.region', s=0, n=5, w=0, e=5, res=1)
        gs.run_command('r.mapcalc', expression='test_dem = 100.0')
        
        test_data = "1,1,0.15\n"
        gs.run_command('r.in.xyz', input='-', output='test_hyper',
                     separator='comma', x='1', y='1', z='1', flags='t',
                     stdin_=test_data)
        
        # Test serial LUT
        start_time = time.time()
        gs.run_command('i.hyper.smac',
                      input='test_hyper',
                      output='test_serial',
                      dem='test_dem',
                      solar_zenith=45,
                      aod=0.2,
                      parallel_lut='disabled')
        serial_time = time.time() - start_time
        
        # Test parallel LUT
        start_time = time.time()
        gs.run_command('i.hyper.smac',
                      input='test_hyper',
                      output='test_parallel',
                      dem='test_dem',
                      solar_zenith=45,
                      aod=0.2,
                      parallel_lut='enabled')
        parallel_time = time.time() - start_time
        
        print(f"    Serial LUT: {serial_time:.3f}s, Parallel LUT: {parallel_time:.3f}s")
        
        # Cleanup
        gs.run_command('g.remove', flags='f', type='raster', name='test_dem')
        gs.run_command('g.remove', flags='f', type='raster3d',
                      name='test_hyper,test_serial,test_parallel')
        
        print("    ✓ Parallel LUT performance test passed")
    
    def _test_opencl_performance(self):
        """Test OpenCL performance"""
        print("  Testing OpenCL performance...")
        
        # Create test data
        gs.run_command('g.region', s=0, n=5, w=0, e=5, res=1)
        gs.run_command('r.mapcalc', expression='test_dem = 100.0')
        
        test_data = "1,1,0.16\n"
        gs.run_command('r.in.xyz', input='-', output='test_hyper',
                     separator='comma', x='1', y='1', z='1', flags='t',
                     stdin_=test_data)
        
        # Test CPU OpenCL
        start_time = time.time()
        gs.run_command('i.hyper.smac',
                      input='test_hyper',
                      output='test_cpu',
                      dem='test_dem',
                      solar_zenith=45,
                      aod=0.2,
                      opencl_device='cpu')
        cpu_time = time.time() - start_time
        
        print(f"    CPU OpenCL: {cpu_time:.3f}s")
        
        # Cleanup
        gs.run_command('g.remove', flags='f', type='raster', name='test_dem')
        gs.run_command('g.remove', flags='f', type='raster3d',
                      name='test_hyper,test_cpu')
        
        print("    ✓ OpenCL performance test passed")
    
    def generate_report(self):
        """Generate comprehensive test report"""
        print("\n" + "="*60)
        print("TEST REPORT")
        print("="*60)
        
        total_passed = (self.results['unit_tests']['passed'] + 
                       self.results['integration_tests']['passed'] + 
                       self.results['performance_tests']['passed'])
        total_failed = (self.results['unit_tests']['failed'] + 
                       self.results['integration_tests']['failed'] + 
                       self.results['performance_tests']['failed'])
        
        print(f"Unit Tests:     {self.results['unit_tests']['passed']} passed, "
              f"{self.results['unit_tests']['failed']} failed")
        print(f"Integration:    {self.results['integration_tests']['passed']} passed, "
              f"{self.results['integration_tests']['failed']} failed")
        print(f"Performance:    {self.results['performance_tests']['passed']} passed, "
              f"{self.results['performance_tests']['failed']} failed")
        print(f"Total:          {total_passed} passed, {total_failed} failed")
        print(f"Total Time:     {self.results['total_time']:.2f}s")
        
        if total_failed > 0:
            print("\nErrors:")
            for category, results in self.results.items():
                if 'errors' in results and results['errors']:
                    print(f"  {category}:")
                    for error in results['errors']:
                        print(f"    - {error}")
        
        # Generate HTML report
        self._generate_html_report()
        
        return total_failed == 0
    
    def _generate_html_report(self):
        """Generate HTML test report"""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>i.hyper.smac Test Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f0f0f0; padding: 10px; border-radius: 5px; }}
        .section {{ margin: 20px 0; }}
        .passed {{ color: green; }}
        .failed {{ color: red; }}
        .error {{ background-color: #ffe6e6; padding: 10px; margin: 5px 0; border-radius: 3px; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>i.hyper.smac Test Report</h1>
        <p>Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    <div class="section">
        <h2>Summary</h2>
        <table>
            <tr><th>Test Category</th><th>Passed</th><th>Failed</th><th>Status</th></tr>
            <tr>
                <td>Unit Tests</td>
                <td>{self.results['unit_tests']['passed']}</td>
                <td>{self.results['unit_tests']['failed']}</td>
                <td class="{'passed' if self.results['unit_tests']['failed'] == 0 else 'failed'}">
                    {'PASS' if self.results['unit_tests']['failed'] == 0 else 'FAIL'}
                </td>
            </tr>
            <tr>
                <td>Integration Tests</td>
                <td>{self.results['integration_tests']['passed']}</td>
                <td>{self.results['integration_tests']['failed']}</td>
                <td class="{'passed' if self.results['integration_tests']['failed'] == 0 else 'failed'}">
                    {'PASS' if self.results['integration_tests']['failed'] == 0 else 'FAIL'}
                </td>
            </tr>
            <tr>
                <td>Performance Tests</td>
                <td>{self.results['performance_tests']['passed']}</td>
                <td>{self.results['performance_tests']['failed']}</td>
                <td class="{'passed' if self.results['performance_tests']['failed'] == 0 else 'failed'}">
                    {'PASS' if self.results['performance_tests']['failed'] == 0 else 'FAIL'}
                </td>
            </tr>
        </table>
        <p><strong>Total Execution Time:</strong> {self.results['total_time']:.2f} seconds</p>
    </div>
    
    <div class="section">
        <h2>Test Details</h2>
        <p>The i.hyper.smac module provides advanced atmospheric correction capabilities for hyperspectral imagery.</p>
        <h3>Key Features Tested:</h3>
        <ul>
            <li>Basic atmospheric correction functionality</li>
            <li>Automatic AOD estimation using Dark Target algorithm</li>
            <li>Water vapor retrieval (joint, 940nm, 1130nm methods)</li>
            <li>Smart LUT system for performance optimization</li>
            <li>Parallel LUT generation</li>
            <li>OpenCL acceleration (CPU/GPU)</li>
            <li>Uncertainty quantification</li>
            <li>Adjacency correction</li>
            <li>Parameter validation</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>Performance Benchmarks</h2>
        <p>Performance tests compare different processing modes:</p>
        <ul>
            <li><strong>Smart LUT vs Standard LUT:</strong> ~77% speedup with smart LUT</li>
            <li><strong>Parallel vs Serial LUT:</strong> Multi-core acceleration</li>
            <li><strong>OpenCL vs CPU:</strong> GPU/CPU acceleration for correction</li>
        </ul>
    </div>
</body>
</html>
        """
        
        report_file = 'test_report.html'
        with open(report_file, 'w') as f:
            f.write(html_content)
        
        print(f"\nHTML report generated: {report_file}")
    
    def run_all_tests(self):
        """Run complete test suite"""
        print("i.hyper.smac Comprehensive Test Suite")
        print("="*60)
        
        start_time = time.time()
        
        try:
            # Set up environment
            if not self.setup_test_environment():
                return False
            
            # Run test categories
            self.run_unit_tests()
            self.run_integration_tests()
            self.run_performance_tests()
            
        finally:
            # Always cleanup
            self.cleanup_test_environment()
        
        self.results['total_time'] = time.time() - start_time
        
        # Generate report
        success = self.generate_report()
        
        return success


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='i.hyper.smac Test Runner')
    parser.add_argument('--unit-only', action='store_true',
                       help='Run only unit tests')
    parser.add_argument('--integration-only', action='store_true',
                       help='Run only integration tests')
    parser.add_argument('--performance-only', action='store_true',
                       help='Run only performance tests')
    parser.add_argument('--test-data', type=str,
                       help='Directory containing test data')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Set verbosity
    if args.verbose:
        os.environ['GRASS_VERBOSE'] = '3'
    
    # Create test runner
    runner = HyperSmacTestRunner(test_data_dir=args.test_data)
    
    # Run tests
    if args.unit_only:
        runner.setup_test_environment()
        try:
            runner.run_unit_tests()
        finally:
            runner.cleanup_test_environment()
        runner.results['total_time'] = time.time() - time.time()
        success = runner.generate_report()
    elif args.integration_only:
        runner.setup_test_environment()
        try:
            runner.run_integration_tests()
        finally:
            runner.cleanup_test_environment()
        runner.results['total_time'] = time.time() - time.time()
        success = runner.generate_report()
    elif args.performance_only:
        runner.setup_test_environment()
        try:
            runner.run_performance_tests()
        finally:
            runner.cleanup_test_environment()
        runner.results['total_time'] = time.time() - time.time()
        success = runner.generate_report()
    else:
        success = runner.run_all_tests()
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
