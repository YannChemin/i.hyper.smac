#!/bin/bash

# Test script to verify the new WVC method functionality

echo "Testing i.hyper.smac with new WVC method parameter..."

# Test 1: Check help output includes new parameter description
echo "=== Test 1: Help output ==="
i.hyper.smac --help | grep -A 3 "water_vapor"

# Test 2: Try using joint method (should work now)
echo ""
echo "=== Test 2: Joint method ==="
echo "Command: i.hyper.smac input=test_rad output=test_out dem=test_dem water_vapor=joint --help"
i.hyper.smac input=test_rad output=test_out dem=test_dem water_vapor=joint --help 2>&1 | head -20

# Test 3: Try using 940nm method
echo ""
echo "=== Test 3: 940nm method ==="
echo "Command: i.hyper.smac input=test_rad output=test_out dem=test_dem water_vapor=940nm --help"
i.hyper.smac input=test_rad output=test_out dem=test_dem water_vapor=940nm --help 2>&1 | head -20

# Test 4: Try using numeric value (original behavior)
echo ""
echo "=== Test 4: Numeric value ==="
echo "Command: i.hyper.smac input=test_rad output=test_out dem=test_dem water_vapor=2.5 --help"
i.hyper.smac input=test_rad output=test_out dem=test_dem water_vapor=2.5 --help 2>&1 | head -20

echo ""
echo "Test completed. Check if water_vapor parameter accepts both methods and numeric values."
