# Create a Python script to extract the band
cat > extract_band.py << 'EOF'
import h5py
import numpy as np
import os

# Input HDF5 file
hdf5_file = "/home/yann/RSDATA/Tanager/Kanpur/20250321_054913_40_4001_basic_radiance_hdf5.h5"
band_number = 110  # 0-based index for band 111

# Output file
output_file = "band_111.raw"

# Open the HDF5 file
with h5py.File(hdf5_file, 'r') as f:
    # Access the dataset
    dataset = f['/HDFEOS/SWATHS/HYP/Data Fields/toa_radiance']
    
    # Extract the specific band (0-based index)
    band_data = dataset[band_number, :, :]
    
    # Save as raw binary file
    band_data.astype('float32').tofile(output_file)
    
    # Also save as text for verification
    np.savetxt('band_111.txt', band_data)
    
    print(f"Band {band_number+1} extracted successfully to {output_file}")
    print(f"Data shape: {band_data.shape}")
    print(f"Min: {np.nanmin(band_data)}, Max: {np.nanmax(band_data)}")
EOF

# Make the script executable
chmod +x extract_band.py

# Run the Python script
python3 extract_band.py

# Now import the band into GRASS
# First, get the region information
g.region -g

# Then import the binary file
r.in.bin input=band_111.raw output=band_111_fixed bytes=4 \
  north=2918190 south=2890950 east=431430 west=406800 \
  rows=732 cols=607 anull=-9999 --o

# Check the imported band
r.info band_111_fixed
r.univar band_111_fixed
