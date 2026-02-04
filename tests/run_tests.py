#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test runner for i.hyper.smac module.

Run with: python tests/run_tests.py

Or using unittest discovery:
    python -m unittest discover -s tests -p 'test_*.py'

Or using pytest:
    pytest tests/
"""

import unittest
import sys
import os

# Add lib directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))


def run_tests():
    """Discover and run all tests."""
    # Get the directory containing this script
    test_dir = os.path.dirname(os.path.abspath(__file__))

    # Create test loader
    loader = unittest.TestLoader()

    # Discover tests
    suite = loader.discover(test_dir, pattern='test_*.py')

    # Create runner with verbosity
    runner = unittest.TextTestRunner(verbosity=2)

    # Run tests
    result = runner.run(suite)

    # Return exit code based on success
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(run_tests())
