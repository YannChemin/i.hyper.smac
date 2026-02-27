# Smart LUT Integration for i.hyper.smac

## Dynamic LUT Coverage Checking and Regeneration

### Problem Solved
Instead of generating a generic LUT with 13 AOD points (156 uvspec runs), the system now:
1. Estimates scene AOD from the image
2. Checks if existing LUT covers the scene AOD adequately  
3. Generates scene-specific LUT only when needed
4. Uses optimal AOD grid (3-5 points instead of 13)

### Performance Improvement
```
Current: 13 AOD × 6 H2O × 2 albedo = 156 uvspec runs (~65 minutes)
Smart:   3 AOD × 6 H2O × 2 albedo = 36 uvspec runs (~15 minutes)
Speedup: 77% reduction in LUT generation time
```

### Integration Code

Add to i.hyper.smac.py in the LUT generation section:

```python
# After AOD estimation (around line 1690)
scene_aod = estimated_aod  # From dark target retrieval

# Smart LUT generation with dynamic coverage checking
try:
    from . import smart_lut
    gs.message("Using smart LUT management with dynamic coverage checking")
    
    # Generate or load optimal LUT for this specific scene
    atm_lut = smart_lut.get_smart_lut_or_generate(
        scene_aod=scene_aod,
        precision_threshold=0.05,  # 5% precision requirement
        sza=solar_zenith, vza=view_zenith, phi=view_azimuth,
        pressure=pressure, aerosol_model=aerosol_model,
        h2o=water_vapor, o3=ozone,
        wl_min=wl_min, wl_max=wl_max, wl_step=2,
        angstrom_alpha=angstrom_alpha,
        parallel=True,
        max_workers=6
    )
    
except ImportError:
    gs.warning("Smart LUT not available, falling back to standard LUT generation")
    # Fallback to existing LUT generation
    atm_lut = lut_module.AtmosphericLUT.get_or_generate(
        sza=solar_zenith, vza=view_zenith, phi=view_azimuth,
        pressure=pressure, aerosol_model=aerosol_model,
        wl_min=wl_min, wl_max=wl_max, wl_step=2,
        h2o=water_vapor, o3=ozone,
        angstrom_alpha=angstrom_alpha,
        force_regenerate=force_regenerate,
        parallel=parallel,
        max_workers=max_workers
    )
```

### Smart LUT Behavior Examples

#### Example 1: Perfect Coverage
```
Scene AOD: 0.664
Existing LUT: [0.5, 0.65, 0.8]
Result: Use existing LUT (no regeneration needed)
```

#### Example 2: Insufficient Precision  
```
Scene AOD: 0.664
Existing LUT: [0.5, 0.8]  # 0.164 away from scene AOD
Precision threshold: 0.05
Result: Regenerate with [0.5, 0.664, 0.8]
```

#### Example 3: Outside Bounds
```
Scene AOD: 1.8
Existing LUT: [0.5, 0.65, 0.8]
Result: Regenerate with [1.26, 1.8, 2.34]
```

### Cache Management

The smart LUT system uses intelligent caching:

```
Cache key: sza40.0_vza50.0_phi30.0_p998.4_urban_0.664_2.00_356-2500_0.05.lut
```

- **Scene-specific**: Different AOD = different cache file
- **Precision-aware**: Different precision = different cache file  
- **Geometry-specific**: Different angles = different cache file
- **Reusable**: Same scene parameters = cached LUT reused

### Configuration Options

```python
# Precision requirements
precision_threshold = 0.05  # 5% AOD precision (default: 0.1 = 10%)

# AOD grid density
n_aod_points = 3  # Minimum bracketing (default: 5 for more precision)
margin_factor = 0.3  # ±30% around scene AOD (default: 0.3)

# Caching
base_cache_dir = "/tmp/grass_lut_cache"  # Custom cache location
```

### Benefits

1. **77% faster LUT generation** for typical scenes
2. **Better precision** around actual scene AOD
3. **Intelligent caching** reduces redundant computations
4. **Automatic fallback** to standard LUT if needed
5. **Dynamic adaptation** to scene characteristics

### Implementation Steps

1. Add `smart_lut.py` to the module list in Makefile
2. Import and integrate smart LUT in i.hyper.smac.py
3. Add configuration options for precision thresholds
4. Test with various AOD ranges and scenes
5. Monitor cache usage and performance

This system makes LUT generation truly dynamic and scene-adaptive!
