"""
OpenCL Accelerator for i.hyper.smac atmospheric correction

This module provides GPU and multi-core acceleration for atmospheric correction
processing using pyOpenCL. It includes device detection, kernel compilation,
and optimized processing for LUT interpolation and atmospheric correction.

Author: Enhanced i.hyper.smac team
"""

import numpy as np
import sys
import os
import time
from typing import Tuple, Optional, Dict, Any, List

# Try to import pyOpenCL, provide fallback if not available
try:
    import pyopencl as cl
    PYOPENCL_AVAILABLE = True
except ImportError:
    PYOPENCL_AVAILABLE = False
    print("Warning: pyOpenCL not available. CPU processing will be used.")
    print("Install with: pip install pyopencl")

# OpenCL kernel source code
OPENCL_KERNELS = """
// Atmospheric correction kernel for hyperspectral imagery
__kernel void atmospheric_correction(
    __global const float* input_data,
    __global float* output_data,
    __global const float* lut_r_atm,
    __global const float* lut_t_down,
    __global const float* lut_t_up,
    __global const float* lut_s,
    __global const float* lut_wavelengths,
    __global const float* lut_aod,
    __global const float* lut_wvc,
    __global const float* lut_o3,
    __global const float* lut_pressure,
    __global const float* lut_sza,
    __global const float* lut_vza,
    __global const float* lut_saa,
    __global const float* lut_vaa,
    __global const int* band_indices,
    __global const float* band_wavelengths,
    __global const float* aod_map,
    __global const float* wvc_map,
    __global const float* o3_map,
    __global const float* pressure_map,
    const float sza,
    const float vza,
    const float saa,
    const float vaa,
    const int num_bands,
    const int num_pixels,
    const int lut_size,
    const int use_maps)
{
    // Get pixel and band indices
    int pixel_id = get_global_id(0);
    int band_id = get_global_id(1);
    
    if (pixel_id >= num_pixels || band_id >= num_bands)
        return;
    
    // Calculate linear index
    int data_idx = pixel_id * num_bands + band_id;
    
    // Get input radiance
    float radiance = input_data[data_idx];
    
    // Get atmospheric parameters for this pixel
    float pixel_aod = use_maps ? aod_map[pixel_id] : aod_map[0];
    float pixel_wvc = use_maps ? wvc_map[pixel_id] : wvc_map[0];
    float pixel_o3 = use_maps ? o3_map[pixel_id] : o3_map[0];
    float pixel_pressure = use_maps ? pressure_map[pixel_id] : pressure_map[0];
    
    // Get band wavelength
    float wavelength = band_wavelengths[band_id];
    
    // Find nearest LUT entries (simplified trilinear interpolation)
    int aod_idx = (int)(pixel_aod * (lut_size - 1) / 2.0f);
    int wvc_idx = (int)(pixel_wvc * (lut_size - 1) / 5.0f);
    int o3_idx = (int)(pixel_o3 * (lut_size - 1) / 1.0f);
    
    // Clamp indices
    aod_idx = clamp(aod_idx, 0, lut_size - 1);
    wvc_idx = clamp(wvc_idx, 0, lut_size - 1);
    o3_idx = clamp(o3_idx, 0, lut_size - 1);
    
    // Calculate LUT index (simplified)
    int lut_idx = aod_idx * lut_size * lut_size + wvc_idx * lut_size + o3_idx;
    
    // Get atmospheric parameters from LUT (simplified - would need proper 6D interpolation)
    float r_atm = lut_r_atm[lut_idx * num_bands + band_id];
    float t_down = lut_t_down[lut_idx * num_bands + band_id];
    float t_up = lut_t_up[lut_idx * num_bands + band_id];
    float s = lut_s[lut_idx * num_bands + band_id];
    
    // Apply atmospheric correction (simplified SMAC equation)
    float denominator = 1.0f + s * radiance;
    if (denominator > 0.001f) {
        output_data[data_idx] = (radiance - r_atm) / (t_down * t_up * denominator);
    } else {
        output_data[data_idx] = 0.0f;
    }
}

// Gas absorption correction kernel
__kernel void gas_absorption_correction(
    __global float* reflectance_data,
    __global const float* gas_transmittance,
    __global const int* band_indices,
    const int num_bands,
    const int num_pixels)
{
    int pixel_id = get_global_id(0);
    int band_id = get_global_id(1);
    
    if (pixel_id >= num_pixels || band_id >= num_bands)
        return;
    
    int data_idx = pixel_id * num_bands + band_id;
    float transmittance = gas_transmittance[band_id];
    
    // Apply gas absorption correction
    if (transmittance > 0.001f) {
        reflectance_data[data_idx] /= transmittance;
    }
}

// Spectral polishing kernel
__kernel void spectral_polishing(
    __global float* reflectance_data,
    __global const float* band_weights,
    const int num_bands,
    const int num_pixels,
    const float threshold)
{
    int pixel_id = get_global_id(0);
    int band_id = get_global_id(1);
    
    if (pixel_id >= num_pixels || band_id >= num_bands)
        return;
    
    int data_idx = pixel_id * num_bands + band_id;
    float value = reflectance_data[data_idx];
    float weight = band_weights[band_id];
    
    // Apply weighted smoothing
    if (weight < threshold) {
        // Simple neighboring band averaging (would need proper implementation)
        float sum = 0.0f;
        float count = 0.0f;
        
        if (band_id > 0) {
            sum += reflectance_data[pixel_id * num_bands + band_id - 1];
            count += 1.0f;
        }
        if (band_id < num_bands - 1) {
            sum += reflectance_data[pixel_id * num_bands + band_id + 1];
            count += 1.0f;
        }
        
        if (count > 0.0f) {
            reflectance_data[data_idx] = sum / count;
        }
    }
}

// Uncertainty calculation kernel
__kernel void calculate_uncertainty(
    __global const float* input_data,
    __global const float* output_data,
    __global float* uncertainty_data,
    __global const float* lut_uncertainty,
    __global const int* band_indices,
    const int num_bands,
    const int num_pixels)
{
    int pixel_id = get_global_id(0);
    int band_id = get_global_id(1);
    
    if (pixel_id >= num_pixels || band_id >= num_bands)
        return;
    
    int data_idx = pixel_id * num_bands + band_id;
    
    // Simplified uncertainty calculation
    float input_uncertainty = 0.01f * input_data[data_idx];  // 1% input uncertainty
    float lut_unc = lut_uncertainty[band_id];  // LUT uncertainty
    float output_uncertainty = sqrt(input_uncertainty * input_uncertainty + lut_unc * lut_unc);
    
    uncertainty_data[data_idx] = output_uncertainty;
}
"""


class OpenCLAccelerator:
    """
    OpenCL accelerator for atmospheric correction processing.
    
    Provides GPU and multi-core acceleration for LUT interpolation,
    atmospheric correction, and post-processing operations.
    """
    
    def __init__(self, device_type='auto', verbose=False):
        """
        Initialize OpenCL accelerator.
        
        Args:
            device_type: 'gpu', 'cpu', 'auto' (default)
            verbose: Enable verbose output
        """
        self.verbose = verbose
        self.ctx = None
        self.queue = None
        self.device = None
        self.program = None
        self.available = PYOPENCL_AVAILABLE
        
        if not PYOPENCL_AVAILABLE:
            return
        
        # Initialize OpenCL context and device
        self._initialize_opencl(device_type)
        
        # Compile kernels
        if self.ctx is not None:
            self._compile_kernels()
    
    def _initialize_opencl(self, device_type):
        """Initialize OpenCL context and select device."""
        try:
            platforms = cl.get_platforms()
            if not platforms:
                if self.verbose:
                    print("No OpenCL platforms found")
                return
            
            # Find suitable device
            devices = []
            for platform in platforms:
                if device_type == 'gpu':
                    platform_devices = platform.get_devices(device_type=cl.device_type.GPU)
                elif device_type == 'cpu':
                    platform_devices = platform.get_devices(device_type=cl.device_type.CPU)
                else:  # auto
                    platform_devices = platform.get_devices()
                devices.extend(platform_devices)
            
            if not devices:
                if self.verbose:
                    print("No suitable OpenCL devices found")
                return
            
            # Select best device (prefer GPU)
            self.device = self._select_best_device(devices)
            if self.device is None:
                return
            
            self.ctx = cl.Context([self.device])
            self.queue = cl.CommandQueue(self.ctx)
            
            if self.verbose:
                print(f"Using OpenCL device: {self.device.name.strip()}")
                print(f"Device type: {self._get_device_type(self.device)}")
                print(f"Max compute units: {self.device.max_compute_units}")
                print(f"Max work group size: {self.device.max_work_group_size}")
                print(f"Global memory: {self.device.global_mem_size // (1024*1024)} MB")
                
        except Exception as e:
            if self.verbose:
                print(f"OpenCL initialization failed: {e}")
            self.available = False
    
    def _select_best_device(self, devices):
        """Select the best available device."""
        # Prefer GPU devices, then CPU
        gpu_devices = [d for d in devices if d.type & cl.device_type.GPU]
        cpu_devices = [d for d in devices if d.type & cl.device_type.CPU]
        
        if gpu_devices:
            # Select GPU with most compute units
            return max(gpu_devices, key=lambda d: d.max_compute_units)
        elif cpu_devices:
            # Select CPU with most compute units
            return max(cpu_devices, key=lambda d: d.max_compute_units)
        else:
            # Return first available device
            return devices[0] if devices else None
    
    def _get_device_type(self, device):
        """Get human-readable device type."""
        if device.type & cl.device_type.GPU:
            return "GPU"
        elif device.type & cl.device_type.CPU:
            return "CPU"
        elif device.type & cl.device_type.ACCELERATOR:
            return "Accelerator"
        else:
            return "Unknown"
    
    def _compile_kernels(self):
        """Compile OpenCL kernels."""
        try:
            self.program = cl.Program(self.ctx, OPENCL_KERNELS).build()
            if self.verbose:
                print("OpenCL kernels compiled successfully")
        except Exception as e:
            if self.verbose:
                print(f"Kernel compilation failed: {e}")
            self.available = False
    
    def is_available(self):
        """Check if OpenCL acceleration is available."""
        return self.available and self.ctx is not None and self.program is not None
    
    def get_device_info(self):
        """Get information about the selected device."""
        if not self.is_available():
            return None
        
        return {
            'name': self.device.name.strip(),
            'type': self._get_device_type(self.device),
            'max_compute_units': self.device.max_compute_units,
            'max_work_group_size': self.device.max_work_group_size,
            'global_memory': self.device.global_mem_size,
            'local_memory': self.device.local_mem_size,
            'extensions': self.device.extensions.strip().split()
        }
    
    def atmospheric_correction_gpu(self, input_data, lut_data, atmospheric_maps, 
                                 geometry, bands_info, use_maps=True):
        """
        Perform atmospheric correction using GPU acceleration.
        
        Args:
            input_data: Input radiance data (pixels x bands)
            lut_data: LUT atmospheric parameters
            atmospheric_maps: AOD, WVC, O3, pressure maps
            geometry: Solar and view angles
            bands_info: Band wavelength information
            use_maps: Whether to use atmospheric maps
            
        Returns:
            Corrected reflectance data
        """
        if not self.is_available():
            return None
        
        try:
            # Handle 3D hyperspectral data (bands, rows, cols)
            if len(input_data.shape) == 3:
                num_bands, num_rows, num_cols = input_data.shape
                num_pixels = num_rows * num_cols
                # Reshape to 2D for GPU processing
                input_data_2d = input_data.reshape(num_bands, -1).T
            else:
                # 2D data (pixels, bands)
                num_pixels, num_bands = input_data.shape
                input_data_2d = input_data
            
            # Prepare data for GPU
            input_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                   hostbuf=input_data_2d.astype(np.float32))
            output_buffer = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY,
                                    size=input_data.nbytes)
            
            # Prepare LUT data (simplified - would need proper LUT structure)
            lut_r_atm = lut_data.get('r_atm', np.zeros((1, num_bands)))
            lut_t_down = lut_data.get('t_down', np.zeros((1, num_bands)))
            lut_t_up = lut_data.get('t_up', np.zeros((1, num_bands)))
            lut_s = lut_data.get('s', np.zeros((1, num_bands)))
            
            # Create LUT buffers
            lut_r_atm_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                        hostbuf=lut_r_atm.astype(np.float32))
            lut_t_down_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                         hostbuf=lut_t_down.astype(np.float32))
            lut_t_up_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                       hostbuf=lut_t_up.astype(np.float32))
            lut_s_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                   hostbuf=lut_s.astype(np.float32))
            
            # Prepare atmospheric maps
            aod_map = atmospheric_maps.get('aod', np.array([0.15]))
            wvc_map = atmospheric_maps.get('wvc', np.array([2.0]))
            o3_map = atmospheric_maps.get('o3', np.array([0.3]))
            pressure_map = atmospheric_maps.get('pressure', np.array([1013.0]))
            
            aod_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                  hostbuf=aod_map.astype(np.float32))
            wvc_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                  hostbuf=wvc_map.astype(np.float32))
            o3_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                hostbuf=o3_map.astype(np.float32))
            pressure_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                      hostbuf=pressure_map.astype(np.float32))
            
            # Band information
            band_wavelengths = np.array([band['wavelength'] for band in bands_info], dtype=np.float32)
            band_wavelengths_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                              hostbuf=band_wavelengths)
            
            # Execute kernel
            global_size = (num_pixels, num_bands)
            local_size = None  # Let OpenCL decide
            
            self.program.atmospheric_correction(
                self.queue, global_size, local_size,
                input_buffer, output_buffer,
                lut_r_atm_buffer, lut_t_down_buffer, lut_t_up_buffer, lut_s_buffer,
                np.zeros(1, dtype=np.float32),  # Dummy lut_wavelengths
                np.zeros(1, dtype=np.float32),  # Dummy lut_aod
                np.zeros(1, dtype=np.float32),  # Dummy lut_wvc
                np.zeros(1, dtype=np.float32),  # Dummy lut_o3
                np.zeros(1, dtype=np.float32),  # Dummy lut_pressure
                np.zeros(1, dtype=np.float32),  # Dummy lut_sza
                np.zeros(1, dtype=np.float32),  # Dummy lut_vza
                np.zeros(1, dtype=np.float32),  # Dummy lut_saa
                np.zeros(1, dtype=np.float32),  # Dummy lut_vaa
                np.arange(num_bands, dtype=np.int32),  # band_indices
                band_wavelengths_buffer,
                aod_buffer, wvc_buffer, o3_buffer, pressure_buffer,
                np.float32(geometry.get('solar_zenith', 40.0)),
                np.float32(geometry.get('view_zenith', 0.0)),
                np.float32(geometry.get('solar_azimuth', 0.0)),
                np.float32(geometry.get('view_azimuth', 0.0)),
                np.int32(num_bands),
                np.int32(num_pixels),
                np.int32(1),  # lut_size (simplified)
                np.int32(1 if use_maps else 0)
            )
            
            # Read results
            output_data_2d = np.empty_like(input_data_2d, dtype=np.float32)
            cl.enqueue_copy(self.queue, output_data_2d, output_buffer)
            
            # Reshape back to original format
            if len(input_data.shape) == 3:
                output_data = output_data_2d.T.reshape(num_bands, num_rows, num_cols)
            else:
                output_data = output_data_2d
            
            return output_data
            
        except Exception as e:
            if self.verbose:
                print(f"GPU atmospheric correction failed: {e}")
            return None
    
    def gas_absorption_correction_gpu(self, reflectance_data, gas_transmittance, bands_info):
        """
        Apply gas absorption correction using GPU acceleration.
        
        Args:
            reflectance_data: Reflectance data (pixels x bands)
            gas_transmittance: Gas transmittance per band
            bands_info: Band wavelength information
            
        Returns:
            Gas-corrected reflectance data
        """
        if not self.is_available():
            return None
        
        try:
            # Handle 3D hyperspectral data (bands, rows, cols)
            if len(reflectance_data.shape) == 3:
                num_bands, num_rows, num_cols = reflectance_data.shape
                num_pixels = num_rows * num_cols
                # Reshape to 2D for GPU processing
                reflectance_data_2d = reflectance_data.reshape(num_bands, -1).T
            else:
                # 2D data (pixels, bands)
                num_pixels, num_bands = reflectance_data.shape
                reflectance_data_2d = reflectance_data
            
            # Prepare buffers
            reflectance_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR,
                                         hostbuf=reflectance_data_2d.astype(np.float32))
            transmittance_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                           hostbuf=gas_transmittance.astype(np.float32))
            
            # Execute kernel
            global_size = (num_pixels, num_bands)
            self.program.gas_absorption_correction(
                self.queue, global_size, None,
                reflectance_buffer, transmittance_buffer,
                np.arange(num_bands, dtype=np.int32),
                np.int32(num_bands),
                np.int32(num_pixels)
            )
            
            # Read results (in-place modification)
            cl.enqueue_copy(self.queue, reflectance_data, reflectance_buffer)
            
            return reflectance_data
            
        except Exception as e:
            if self.verbose:
                print(f"GPU gas absorption correction failed: {e}")
            return None
    
    def spectral_polishing_gpu(self, reflectance_data, band_weights, threshold=0.1):
        """
        Apply spectral polishing using GPU acceleration.
        
        Args:
            reflectance_data: Reflectance data (pixels x bands)
            band_weights: Quality weights per band
            threshold: Weight threshold for polishing
            
        Returns:
            Polished reflectance data
        """
        if not self.is_available():
            return None
        
        try:
            # Handle 3D hyperspectral data (bands, rows, cols)
            if len(reflectance_data.shape) == 3:
                num_bands, num_rows, num_cols = reflectance_data.shape
                num_pixels = num_rows * num_cols
                # Reshape to 2D for GPU processing
                reflectance_data_2d = reflectance_data.reshape(num_bands, -1).T
            else:
                # 2D data (pixels, bands)
                num_pixels, num_bands = reflectance_data.shape
                reflectance_data_2d = reflectance_data
            
            # Prepare buffers
            reflectance_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR,
                                         hostbuf=reflectance_data_2d.astype(np.float32))
            weights_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                      hostbuf=band_weights.astype(np.float32))
            
            # Execute kernel
            global_size = (num_pixels, num_bands)
            self.program.spectral_polishing(
                self.queue, global_size, None,
                reflectance_buffer, weights_buffer,
                np.int32(num_bands),
                np.int32(num_pixels),
                np.float32(threshold)
            )
            
            # Read results (in-place modification)
            cl.enqueue_copy(self.queue, reflectance_data, reflectance_buffer)
            
            return reflectance_data
            
        except Exception as e:
            if self.verbose:
                print(f"GPU spectral polishing failed: {e}")
            return None
    
    def calculate_uncertainty_gpu(self, input_data, output_data, lut_uncertainty, bands_info):
        """
        Calculate uncertainty using GPU acceleration.
        
        Args:
            input_data: Input radiance data
            output_data: Output reflectance data
            lut_uncertainty: LUT uncertainty per band
            bands_info: Band wavelength information
            
        Returns:
            Uncertainty data
        """
        if not self.is_available():
            return None
        
        try:
            # Handle 3D hyperspectral data (bands, rows, cols)
            if len(input_data.shape) == 3:
                num_bands, num_rows, num_cols = input_data.shape
                num_pixels = num_rows * num_cols
                # Reshape to 2D for GPU processing
                input_data_2d = input_data.reshape(num_bands, -1).T
            else:
                # 2D data (pixels, bands)
                num_pixels, num_bands = input_data.shape
                input_data_2d = input_data
            
            # Prepare buffers
            input_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                   hostbuf=input_data_2d.astype(np.float32))
            output_buffer = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY,
                                    size=input_data.nbytes)
            uncertainty_buffer = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY,
                                          size=input_data.nbytes)
            lut_unc_buffer = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                     hostbuf=lut_uncertainty.astype(np.float32))
            
            # Execute kernel
            global_size = (num_pixels, num_bands)
            self.program.calculate_uncertainty(
                self.queue, global_size, None,
                input_buffer, output_buffer, uncertainty_buffer, lut_unc_buffer,
                np.arange(num_bands, dtype=np.int32),
                np.int32(num_bands),
                np.int32(num_pixels)
            )
            
            # Read results
            uncertainty_data = np.empty_like(input_data, dtype=np.float32)
            cl.enqueue_copy(self.queue, uncertainty_data, uncertainty_buffer)
            
            return uncertainty_data
            
        except Exception as e:
            if self.verbose:
                print(f"GPU uncertainty calculation failed: {e}")
            return None
    
    def benchmark_performance(self, test_data_size=(1000, 100)):
        """
        Benchmark OpenCL performance vs CPU.
        
        Args:
            test_data_size: Tuple of (num_pixels, num_bands) for testing
            
        Returns:
            Dictionary with performance metrics
        """
        if not self.is_available():
            return {'error': 'OpenCL not available'}
        
        try:
            num_pixels, num_bands = test_data_size
            
            # Generate test data
            input_data = np.random.rand(num_pixels, num_bands).astype(np.float32)
            gas_transmittance = np.random.rand(num_bands).astype(np.float32)
            band_weights = np.random.rand(num_bands).astype(np.float32)
            
            # Test GPU performance
            start_time = time.time()
            gpu_result = self.gas_absorption_correction_gpu(
                input_data.copy(), gas_transmittance, [])
            gpu_time = time.time() - start_time
            
            # Test CPU performance (simple loop)
            start_time = time.time()
            cpu_result = input_data.copy() / gas_transmittance[np.newaxis, :]
            cpu_time = time.time() - start_time
            
            # Calculate speedup
            speedup = cpu_time / gpu_time if gpu_time > 0 else 0
            
            return {
                'device_info': self.get_device_info(),
                'test_size': test_data_size,
                'gpu_time': gpu_time,
                'cpu_time': cpu_time,
                'speedup': speedup,
                'gpu_available': True
            }
            
        except Exception as e:
            return {'error': str(e), 'gpu_available': False}
    
    def cleanup(self):
        """Clean up OpenCL resources."""
        if self.ctx is not None:
            # OpenCL resources are automatically cleaned up when context is garbage collected
            self.ctx = None
            self.queue = None
            self.device = None
            self.program = None


def create_opencl_accelerator(device_type='auto', verbose=False):
    """
    Create and return an OpenCL accelerator instance.
    
    Args:
        device_type: 'gpu', 'cpu', 'auto'
        verbose: Enable verbose output
        
    Returns:
        OpenCLAccelerator instance
    """
    return OpenCLAccelerator(device_type=device_type, verbose=verbose)


def detect_opencl_devices():
    """
    Detect available OpenCL devices.
    
    Returns:
        List of device information dictionaries
    """
    if not PYOPENCL_AVAILABLE:
        return []
    
    try:
        platforms = cl.get_platforms()
        devices_info = []
        
        for platform in platforms:
            devices = platform.get_devices()
            for device in devices:
                device_info = {
                    'platform': platform.name.strip(),
                    'name': device.name.strip(),
                    'type': 'GPU' if device.type & cl.device_type.GPU else 
                           'CPU' if device.type & cl.device_type.CPU else 'Other',
                    'max_compute_units': device.max_compute_units,
                    'max_work_group_size': device.max_work_group_size,
                    'global_memory': device.global_mem_size,
                    'local_memory': device.local_mem_size,
                    'available': device.available
                }
                devices_info.append(device_info)
        
        return devices_info
        
    except Exception as e:
        print(f"Error detecting OpenCL devices: {e}")
        return []


if __name__ == "__main__":
    # Test OpenCL accelerator
    print("Testing OpenCL Accelerator...")
    
    # Detect devices
    devices = detect_opencl_devices()
    print(f"Found {len(devices)} OpenCL devices:")
    for i, device in enumerate(devices):
        print(f"  {i+1}. {device['name']} ({device['type']}) - {device['max_compute_units']} compute units")
    
    # Create accelerator
    accelerator = create_opencl_accelerator(verbose=True)
    
    if accelerator.is_available():
        print("OpenCL accelerator is available!")
        
        # Benchmark performance
        benchmark = accelerator.benchmark_performance()
        if 'error' not in benchmark:
            print(f"Performance test: GPU {benchmark['gpu_time']:.3f}s vs CPU {benchmark['cpu_time']:.3f}s")
            print(f"Speedup: {benchmark['speedup']:.2f}x")
        else:
            print(f"Benchmark failed: {benchmark['error']}")
    else:
        print("OpenCL accelerator not available - CPU processing will be used")
    
    accelerator.cleanup()
