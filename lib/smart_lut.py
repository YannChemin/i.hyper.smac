"""
Smart LUT management with dynamic coverage checking and regeneration.

This module provides intelligent LUT management that:
1. Checks if existing LUT covers scene AOD range adequately
2. Calculates optimal AOD grid for scene-specific processing
3. Triggers regeneration when coverage is insufficient
4. Maintains cache efficiency with smart naming
"""

import numpy as np
import os
import grass.script as gs
from typing import Tuple, Optional, List
import hashlib


class SmartLUTManager:
    """
    Intelligent LUT management with dynamic coverage checking.
    """
    
    def __init__(self, base_cache_dir: str = "/tmp/grass_lut_cache"):
        """
        Initialize smart LUT manager.
        
        Args:
            base_cache_dir: Base directory for LUT caching
        """
        self.base_cache_dir = base_cache_dir
        os.makedirs(base_cache_dir, exist_ok=True)
    
    def check_aod_coverage(self, scene_aod: float, lut_aod_grid: np.ndarray, 
                         precision_threshold: float = 0.1) -> Tuple[bool, str]:
        """
        Check if LUT adequately covers the scene AOD.
        
        Args:
            scene_aod: Estimated AOD from scene
            lut_aod_grid: AOD values in existing LUT
            precision_threshold: Maximum allowed distance to nearest LUT point
            
        Returns:
            Tuple of (coverage_ok, reason)
        """
        # Check if scene AOD is within LUT bounds
        if scene_aod < lut_aod_grid.min():
            return False, f"Scene AOD {scene_aod:.3f} below LUT minimum {lut_aod_grid.min():.3f}"
        
        if scene_aod > lut_aod_grid.max():
            return False, f"Scene AOD {scene_aod:.3f} above LUT maximum {lut_aod_grid.max():.3f}"
        
        # Check precision: distance to nearest LUT point
        distances = np.abs(lut_aod_grid - scene_aod)
        min_distance = distances.min()
        
        if min_distance > precision_threshold:
            return False, f"Scene AOD {scene_aod:.3f} too far from nearest LUT point ({min_distance:.3f} > {precision_threshold})"
        
        return True, f"AOD coverage adequate (nearest point: {min_distance:.3f})"
    
    def calculate_optimal_aod_grid(self, scene_aod: float, 
                                 n_points: int = 5,
                                 margin_factor: float = 0.3) -> np.ndarray:
        """
        Calculate optimal AOD grid for scene-specific processing.
        
        Args:
            scene_aod: Estimated AOD from scene
            n_points: Number of AOD points to generate
            margin_factor: Margin factor around scene AOD (0.3 = ±30%)
            
        Returns:
            Optimal AOD grid array
        """
        # Calculate range around scene AOD
        aod_range = scene_aod * margin_factor
        aod_min = max(0.001, scene_aod - aod_range)  # Minimum AOD > 0
        aod_max = scene_aod + aod_range
        
        # Generate evenly spaced points
        if n_points == 3:
            # Simple bracketing: lower, scene, upper
            return np.array([aod_min, scene_aod, aod_max])
        else:
            # Even spacing with scene AOD included
            grid = np.linspace(aod_min, aod_max, n_points)
            
            # Ensure scene AOD is included (replace nearest point)
            distances = np.abs(grid - scene_aod)
            scene_idx = distances.argmin()
            grid[scene_idx] = scene_aod
            
            return np.sort(grid)
    
    def generate_cache_key(self, sza: float, vza: float, phi: float,
                          pressure: float, aerosol_model: str,
                          scene_aod: float, h2o: float,
                          wl_range: str, precision_threshold: float) -> str:
        """
        Generate cache key for LUT with coverage parameters.
        
        Returns:
            Cache filename hash
        """
        key_data = f"{sza:.1f}_{vza:.1f}_{phi:.1f}_{pressure:.1f}_{aerosol_model}_{scene_aod:.3f}_{h2o:.2f}_{wl_range}_{precision_threshold}"
        return hashlib.md5(key_data.encode()).hexdigest() + ".lut"
    
    def get_smart_lut_config(self, scene_aod: float, precision_threshold: float = 0.1) -> dict:
        """
        Get smart LUT configuration for scene.
        
        Args:
            scene_aod: Estimated AOD from scene
            precision_threshold: Required precision for AOD coverage
            
        Returns:
            Dictionary with LUT configuration
        """
        # Calculate optimal AOD grid
        optimal_aod_grid = self.calculate_optimal_aod_grid(scene_aod)
        
        return {
            'aod_grid': optimal_aod_grid,
            'scene_aod': scene_aod,
            'precision_threshold': precision_threshold,
            'n_aod_points': len(optimal_aod_grid),
            'coverage_range': (optimal_aod_grid.min(), optimal_aod_grid.max())
        }
    
    def should_regenerate_lut(self, scene_aod: float, existing_lut_path: str,
                            precision_threshold: float = 0.1) -> Tuple[bool, str]:
        """
        Determine if LUT should be regenerated.
        
        Args:
            scene_aod: Estimated AOD from scene
            existing_lut_path: Path to existing LUT file
            precision_threshold: Required precision for AOD coverage
            
        Returns:
            Tuple of (should_regenerate, reason)
        """
        if not os.path.exists(existing_lut_path):
            return True, "No existing LUT found"
        
        try:
            # Load existing LUT and check coverage
            from . import lut
            existing_lut = lut.AtmosphericLUT.load(existing_lut_path)
            
            coverage_ok, reason = self.check_aod_coverage(
                scene_aod, existing_lut.aod_grid, precision_threshold
            )
            
            return not coverage_ok, reason
            
        except Exception as e:
            return True, f"Error loading existing LUT: {e}"
    
    def estimate_lut_generation_time(self, n_aod_points: int, n_h2o_points: int,
                                    n_albedo_points: int = 2, 
                                    avg_uvspec_time: float = 25.0) -> dict:
        """
        Estimate LUT generation time and resource requirements.
        
        Args:
            n_aod_points: Number of AOD points in grid
            n_h2o_points: Number of H2O points in grid  
            n_albedo_points: Number of albedo points (usually 2)
            avg_uvspec_time: Average time per uvspec run (seconds)
            
        Returns:
            Dictionary with time estimates
        """
        total_runs = n_aod_points * n_h2o_points * n_albedo_points
        sequential_time = total_runs * avg_uvspec_time
        
        # Parallel speedup estimate (75% resource usage)
        import multiprocessing as mp
        n_workers = max(1, int(mp.cpu_count() * 0.75))
        parallel_time = sequential_time / n_workers
        
        return {
            'total_runs': total_runs,
            'sequential_time_minutes': sequential_time / 60,
            'parallel_time_minutes': parallel_time / 60,
            'n_workers': n_workers,
            'speedup_factor': n_workers
        }


def get_smart_lut_or_generate(scene_aod: float, precision_threshold: float = 0.1,
                              force_regenerate: bool = False, **lut_kwargs) -> 'AtmosphericLUT':
    """
    Get smart LUT with automatic coverage checking and regeneration.
    
    Args:
        scene_aod: Estimated AOD from scene
        precision_threshold: Required precision for AOD coverage
        force_regenerate: Force regeneration regardless of coverage
        **lut_kwargs: Additional arguments for LUT generation
        
    Returns:
        AtmosphericLUT instance with optimal coverage
    """
    manager = SmartLUTManager()
    
    # Get smart LUT configuration
    smart_config = manager.get_smart_lut_config(scene_aod, precision_threshold)
    
    # Generate cache key
    # Extract wavelength parameters and convert to range string
    wl_min = lut_kwargs.get('wl_min', 400)
    wl_max = lut_kwargs.get('wl_max', 2500)
    wl_step = lut_kwargs.get('wl_step', 2)
    wl_range = f"{wl_min}-{wl_max}-{wl_step}"
    
    cache_key = manager.generate_cache_key(
        sza=lut_kwargs.get('sza'),
        vza=lut_kwargs.get('vza'), 
        phi=lut_kwargs.get('phi'),
        pressure=lut_kwargs.get('pressure'),
        aerosol_model=lut_kwargs.get('aerosol_model'),
        scene_aod=scene_aod,
        h2o=lut_kwargs.get('h2o'),
        wl_range=wl_range,
        precision_threshold=precision_threshold
    )
    cache_path = os.path.join(manager.base_cache_dir, cache_key)
    
    # Check if regeneration is needed
    if not force_regenerate:
        should_regenerate, reason = manager.should_regenerate_lut(
            scene_aod, cache_path, precision_threshold
        )
        
        if not should_regenerate:
            gs.message(f"Using existing smart LUT: {cache_path}")
            gs.message(f"LUT coverage: {smart_config['coverage_range'][0]:.3f} - {smart_config['coverage_range'][1]:.3f}")
            from . import lut
            return lut.AtmosphericLUT.load(cache_path)
        else:
            gs.message(f"Regenerating LUT: {reason}")
    else:
        gs.message("Force regenerating smart LUT")
    
    # Generate time estimates
    time_estimates = manager.estimate_lut_generation_time(
        smart_config['n_aod_points'], 
        len(lut_kwargs.get('h2o_grid', [3.5, 5.0, 7.0, 9.5, 12.0, 17.5]))
    )
    
    gs.message(f"Smart LUT generation estimates:")
    gs.message(f"  Total runs: {time_estimates['total_runs']}")
    gs.message(f"  Parallel time: {time_estimates['parallel_time_minutes']:.1f} minutes")
    gs.message(f"  Using {time_estimates['n_workers']} workers")
    
    # Generate new LUT with optimal AOD grid
    import lut
    new_lut = lut.AtmosphericLUT.generate(
        aod_grid=smart_config['aod_grid'],
        **lut_kwargs
    )
    
    # Save with smart cache key
    new_lut.save(cache_path)
    
    gs.message(f"Generated smart LUT: {cache_path}")
    gs.message(f"AOD coverage: {smart_config['coverage_range'][0]:.3f} - {smart_config['coverage_range'][1]:.3f} (±{precision_threshold:.2f})")
    
    return new_lut


# Example usage in i.hyper.smac:
"""
# Replace existing LUT generation with smart version:
scene_aod = estimated_aod  # From dark target retrieval
lut = get_smart_lut_or_generate(
    scene_aod=scene_aod,
    precision_threshold=0.05,  # Tighter precision
    sza=solar_zenith, vza=view_zenith, phi=view_azimuth,
    pressure=pressure, aerosol_model=aerosol_model,
    h2o=water_vapor, o3=ozone,
    parallel=True, max_workers=6
)
"""
