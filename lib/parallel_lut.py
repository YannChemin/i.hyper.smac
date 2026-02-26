"""
Parallel LUT Generation for i.hyper.smac

This module implements parallel uvspec execution to dramatically speed up
LUT generation by running multiple libRadtran instances simultaneously.

Author: Enhanced i.hyper.smac team
"""

import os
import sys
import subprocess
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import tempfile
import numpy as np
import time
from typing import List, Dict, Tuple, Optional
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class ParallelLUTGenerator:
    """
    Parallel LUT generator using multiple uvspec processes.
    
    Dramatically reduces LUT generation time by distributing
    libRadtran runs across available CPU cores.
    """
    
    def __init__(self, max_workers: Optional[int] = None, 
                 uvspec_path: str = "uvspec"):
        """
        Initialize parallel LUT generator.
        
        Args:
            max_workers: Maximum number of parallel processes (default: CPU cores - 1)
            uvspec_path: Path to uvspec executable
        """
        self.max_workers = max_workers or (mp.cpu_count() - 1)
        self.uvspec_path = uvspec_path
        logger.info(f"Initialized parallel LUT generator with {self.max_workers} workers")
    
    def _create_uvspec_input(self, aod: float, h2o: float, albedo: float,
                          sza: float, vza: float, phi: float,
                          pressure: float, aerosol_model: str,
                          output_file: str, wavelengths: str) -> str:
        """
        Create uvspec input file content.
        
        Args:
            aod: Aerosol optical depth at 550nm
            h2o: Water vapor column in g/cm²
            albedo: Surface albedo (0 or 0.5)
            sza: Solar zenith angle in degrees
            vza: View zenith angle in degrees
            phi: Relative azimuth angle in degrees
            pressure: Atmospheric pressure in hPa
            aerosol_model: Aerosol model type
            output_file: Output file path
            wavelengths: Wavelength specification
            
        Returns:
            uvspec input file content as string
        """
        # Convert wavelengths to uvspec format
        if wavelengths.startswith("wl "):
            wl_spec = wavelengths
        else:
            # Convert range like "356-2500" to uvspec format
            if "-" in wavelengths:
                start, end = wavelengths.split("-")
                wl_spec = f"wl {start} {end}"
            else:
                wl_spec = f"wl {wavelengths}"
        
        input_content = f"""\
# Parallel LUT Generation Input
# AOD: {aod:.3f}, H2O: {h2o:.2f}, Albedo: {albedo:.1f}
# Geometry: SZA={sza:.1f}, VZA={vza:.1f}, PHI={phi:.1f}
# Pressure: {pressure:.1f} hPa, Aerosol: {aerosol_model}

atmosphere_file
data_source_file
mol_abs_file
rayleigh_file
aerosol_file
cloud_file 0

wavelength {wl_spec}

output_file {output_file}
radiances

# Advanced options for speed
disort 16
quiet
"""
        return input_content
    
    def _run_single_uvspec(self, params: Tuple) -> Tuple[bool, str, float]:
        """
        Run a single uvspec instance.
        
        Args:
            params: Tuple of (aod, h2o, albedo, sza, vza, phi, pressure, 
                           aerosol_model, output_file, wavelengths)
            
        Returns:
            Tuple of (success, output_file, execution_time)
        """
        (aod, h2o, albedo, sza, vza, phi, pressure, 
         aerosol_model, output_file, wavelengths) = params
        
        start_time = time.time()
        
        try:
            # Create input file
            input_content = self._create_uvspec_input(
                aod, h2o, albedo, sza, vza, phi, pressure,
                aerosol_model, output_file, wavelengths
            )
            
            # Write input file
            input_file = output_file.replace('.out', '.inp')
            with open(input_file, 'w') as f:
                f.write(input_content)
            
            # Run uvspec
            cmd = [self.uvspec_path, input_file]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout per run
            )
            
            execution_time = time.time() - start_time
            
            success = result.returncode == 0
            
            if success:
                logger.info(f"Completed: AOD={aod:.3f}, H2O={h2o:.2f}, "
                           f"Albedo={albedo:.1f} in {execution_time:.1f}s")
            else:
                logger.error(f"Failed: AOD={aod:.3f}, H2O={h2o:.2f}, "
                            f"Albedo={albedo:.1f} - {result.stderr}")
            
            # Cleanup input file
            try:
                os.remove(input_file)
            except:
                pass
            
            return success, output_file, execution_time
            
        except subprocess.TimeoutExpired:
            logger.error(f"Timeout: AOD={aod:.3f}, H2O={h2o:.2f}, Albedo={albedo:.1f}")
            return False, output_file, 300.0
        except Exception as e:
            logger.error(f"Error: AOD={aod:.3f}, H2O={h2o:.2f}, Albedo={albedo:.1f} - {e}")
            return False, output_file, 0.0
    
    def generate_lut_parallel(self, lut_config: Dict) -> bool:
        """
        Generate LUT using parallel uvspec execution.
        
        Args:
            lut_config: Dictionary containing LUT generation parameters:
                - sza: Solar zenith angle
                - vza: View zenith angle  
                - phi: Relative azimuth angle
                - pressure: Atmospheric pressure
                - aerosol_model: Aerosol model
                - aod_values: List of AOD values
                - h2o_values: List of H2O values
                - albedo_values: List of albedo values (0, 0.5)
                - wavelength_range: Wavelength range (e.g., "356-2500")
                - output_prefix: Output file prefix
                - temp_dir: Temporary directory
                
        Returns:
            True if successful, False otherwise
        """
        logger.info("Starting parallel LUT generation...")
        logger.info(f"Configuration: {lut_config}")
        
        # Extract parameters
        sza = lut_config['sza']
        vza = lut_config['vza']
        phi = lut_config.get('phi', 30.0)
        pressure = lut_config['pressure']
        aerosol_model = lut_config['aerosol_model']
        aod_values = lut_config['aod_values']
        h2o_values = lut_config['h2o_values']
        albedo_values = lut_config['albedo_values']
        wavelength_range = lut_config['wavelength_range']
        output_prefix = lut_config['output_prefix']
        temp_dir = lut_config.get('temp_dir', '/tmp')
        
        # Create parameter combinations
        param_combinations = []
        for aod in aod_values:
            for h2o in h2o_values:
                for albedo in albedo_values:
                    output_file = os.path.join(temp_dir, f"{output_prefix}_aod{aod:.2f}_h2o{h2o:.2f}_alb{albedo:.1f}.out")
                    param_combinations.append((aod, h2o, albedo, sza, vza, phi, 
                                           pressure, aerosol_model, output_file, wavelength_range))
        
        total_runs = len(param_combinations)
        logger.info(f"Total LUT runs: {total_runs}")
        logger.info(f"Parallel workers: {self.max_workers}")
        
        # Track progress
        successful_runs = 0
        failed_runs = 0
        total_time = 0.0
        
        # Execute runs in parallel
        start_time = time.time()
        
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all jobs
            future_to_params = {
                executor.submit(self._run_single_uvspec, params): params 
                for params in param_combinations
            }
            
            # Process completed jobs
            for future in as_completed(future_to_params):
                try:
                    success, output_file, execution_time = future.result()
                    total_time += execution_time
                    
                    if success:
                        successful_runs += 1
                    else:
                        failed_runs += 1
                        # Cleanup failed output file
                        try:
                            if os.path.exists(output_file):
                                os.remove(output_file)
                        except:
                            pass
                    
                    # Progress update
                    progress = (successful_runs + failed_runs) / total_runs * 100
                    logger.info(f"Progress: {progress:.1f}% ({successful_runs + failed_runs}/{total_runs}) "
                               f"- Success: {successful_runs}, Failed: {failed_runs}")
                    
                except Exception as e:
                    logger.error(f"Job processing error: {e}")
                    failed_runs += 1
        
        total_execution_time = time.time() - start_time
        
        logger.info(f"Parallel LUT generation completed in {total_execution_time:.1f}s")
        logger.info(f"Total uvspec time: {total_time:.1f}s")
        logger.info(f"Successful runs: {successful_runs}/{total_runs}")
        logger.info(f"Failed runs: {failed_runs}/{total_runs}")
        
        if failed_runs == 0:
            logger.info("✅ All LUT runs completed successfully!")
            return True
        else:
            logger.warning(f"⚠️  {failed_runs} LUT runs failed")
            return False
    
    def estimate_optimal_workers(self, total_runs: int, memory_per_run_mb: float = 100) -> int:
        """
        Estimate optimal number of workers based on available memory and runs.
        
        Args:
            total_runs: Total number of uvspec runs needed
            memory_per_run_mb: Estimated memory usage per run in MB
            
        Returns:
            Optimal number of workers
        """
        # Get available memory (simplified)
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemAvailable:'):
                        available_mb = int(line.split()[1]) // 1024
                        break
                else:
                    # Fallback to total memory
                    for line in f:
                        if line.startswith('MemTotal:'):
                            available_mb = int(line.split()[1]) // 1024 // 4  # Use 25% of total
                            break
        except:
            available_mb = 4096  # 4GB fallback
        
        # Calculate optimal workers
        max_workers_by_memory = int(available_mb / memory_per_run_mb)
        cpu_cores = mp.cpu_count()
        
        optimal_workers = min(
            max_workers_by_memory,
            cpu_cores - 1,  # Leave one core free
            total_runs,
            self.max_workers
        )
        
        logger.info(f"Memory optimization: {available_mb}MB available, "
                   f"{memory_per_run_mb}MB per run → {max_workers_by_memory} workers by memory")
        logger.info(f"CPU optimization: {cpu_cores} cores → {cpu_cores - 1} workers by CPU")
        logger.info(f"Selected optimal workers: {optimal_workers}")
        
        return optimal_workers


def create_parallel_lut_config(sza: float, vza: float, phi: float,
                           pressure: float, aerosol_model: str,
                           wavelength_range: str = "356-2500",
                           output_prefix: str = "lut") -> Dict:
    """
    Create standard LUT configuration for parallel generation.
    
    Args:
        sza: Solar zenith angle in degrees
        vza: View zenith angle in degrees
        phi: Relative azimuth angle in degrees
        pressure: Atmospheric pressure in hPa
        aerosol_model: Aerosol model type
        wavelength_range: Wavelength range specification
        output_prefix: Output file prefix
        
    Returns:
        Dictionary with LUT generation parameters
    """
    return {
        'sza': sza,
        'vza': vza,
        'phi': phi,
        'pressure': pressure,
        'aerosol_model': aerosol_model,
        'aod_values': [0.0, 0.1, 0.3, 0.6, 0.9, 1.2, 1.8, 2.7],
        'h2o_values': [0.3, 0.7, 1.5, 2.5, 3.5, 5.0],
        'albedo_values': [0.0, 0.5],
        'wavelength_range': wavelength_range,
        'output_prefix': output_prefix,
        'temp_dir': tempfile.gettempdir()
    }


def benchmark_sequential_vs_parallel(lut_config: Dict, max_workers: int = 4) -> Dict:
    """
    Benchmark sequential vs parallel LUT generation.
    
    Args:
        lut_config: LUT configuration
        max_workers: Maximum number of parallel workers
        
    Returns:
        Dictionary with benchmark results
    """
    logger.info("Running benchmark...")
    
    # Count total runs
    aod_values = lut_config['aod_values']
    h2o_values = lut_config['h2o_values']
    albedo_values = lut_config['albedo_values']
    total_runs = len(aod_values) * len(h2o_values) * len(albedo_values)
    
    # Sequential timing (estimate)
    estimated_sequential_time = total_runs * 15  # 15 seconds per run estimate
    
    # Parallel timing estimate
    estimated_parallel_time = (total_runs / max_workers) * 15
    
    speedup = estimated_sequential_time / estimated_parallel_time
    
    results = {
        'total_runs': total_runs,
        'sequential_time_estimate': estimated_sequential_time,
        'parallel_time_estimate': estimated_parallel_time,
        'speedup': speedup,
        'max_workers': max_workers,
        'time_saved': estimated_sequential_time - estimated_parallel_time
    }
    
    logger.info(f"Benchmark Results:")
    logger.info(f"  Total runs: {total_runs}")
    logger.info(f"  Sequential estimate: {estimated_sequential_time:.1f}s")
    logger.info(f"  Parallel estimate: {estimated_parallel_time:.1f}s")
    logger.info(f"  Speedup: {speedup:.1f}x")
    logger.info(f"  Time saved: {estimated_parallel_time - estimated_parallel_time:.1f}s")
    
    return results


if __name__ == "__main__":
    # Example usage
    print("Parallel LUT Generator for i.hyper.smac")
    print("=" * 50)
    
    # Example configuration
    config = create_parallel_lut_config(
        sza=40.0, vza=50.0, phi=30.0,
        pressure=1013.25, aerosol_model="urban",
        wavelength_range="356-2500",
        output_prefix="test_lut"
    )
    
    # Create generator
    generator = ParallelLUTGenerator(max_workers=4)
    
    # Run benchmark
    benchmark = benchmark_sequential_vs_parallel(config, max_workers=4)
    
    # Generate LUT
    success = generator.generate_lut_parallel(config)
    
    if success:
        print("✅ Parallel LUT generation completed successfully!")
    else:
        print("❌ Some LUT runs failed")
    
    print(f"Benchmark: {benchmark['speedup']:.1f}x speedup potential")
