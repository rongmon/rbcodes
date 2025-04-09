#!/usr/bin/env python3
"""
Script to restructure the rbcodes repository into a more standard Python package layout.
Run this script from the root of the repository.
"""

import os
import shutil
from pathlib import Path

def create_directory(path):
    """Create directory if it doesn't exist."""
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Created directory: {path}")

def copy_item(source, dest):
    """Copy file or directory, handling existing destinations."""
    # Convert to Path objects if they're not already
    source = Path(source)
    dest = Path(dest)
    
    # Skip if source and destination are the same file
    if source == dest or os.path.samefile(source, dest) if os.path.exists(dest) and os.path.exists(source) else False:
        print(f"Skipping same file: {source}")
        return
        
    if not os.path.exists(source):
        print(f"Warning: Source does not exist: {source}")
        return
        
    if os.path.isdir(source):
        if os.path.exists(dest):
            # If destination exists, copy contents instead of the directory itself
            for item in os.listdir(source):
                s_item = source / item
                d_item = dest / item
                copy_item(s_item, d_item)
            print(f"Merged directory contents: {source} -> {dest}")
        else:
            shutil.copytree(source, dest)
            print(f"Copied directory: {source} -> {dest}")
    else:
        # For files, just overwrite if they exist
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        try:
            shutil.copy2(source, dest)
            print(f"Copied file: {source} -> {dest}")
        except shutil.SameFileError:
            print(f"Skipping same file: {source}")

def main():
    # Define the root and target directories
    root_dir = Path(".")
    src_dir = root_dir / "src" / "rbcodes"
    
    # Create the new directory structure
    create_directory(src_dir)
    create_directory(root_dir / "tests")
    create_directory(root_dir / "docs")
    create_directory(root_dir / "scripts")
    
    # Create subdirectories in src/rbcodes
    subdirs = ["gui", "igm", "catalog", "halo", "lensing", "stats", "utils"]
    for subdir in subdirs:
        create_directory(src_dir / subdir)
    
    # Create test subdirectories
    test_subdirs = [f"test_{subdir}" for subdir in subdirs]
    for test_subdir in test_subdirs:
        create_directory(root_dir / "tests" / test_subdir)
    
    # Create doc subdirectories
    for subdir in subdirs:
        create_directory(root_dir / "docs" / subdir)
    
    # Move files from GUIs to src/rbcodes/gui
    print("Moving GUI files...")
    # Create GUI subdirectories
    gui_subdirs = ["abstools", "multispecviewer", "spectral_gui", "zgui"]
    for subdir in gui_subdirs:
        create_directory(src_dir / "gui" / subdir)
        
        # Move the contents
        source_dir = root_dir / "GUIs" / subdir
        if os.path.exists(source_dir):
            for item in os.listdir(source_dir):
                s = source_dir / item
                d = src_dir / "gui" / subdir / item
                copy_item(s, d)
    
    # Move other GUI files to gui directory
    source_dir = root_dir / "GUIs"
    if os.path.exists(source_dir):
        for item in os.listdir(source_dir):
            if item not in gui_subdirs:  # Skip already processed subdirectories
                s = source_dir / item
                d = src_dir / "gui" / item
                copy_item(s, d)
    
    # Move files from IGM to src/rbcodes/igm
    print("Moving IGM files...")
    igm_dir = root_dir / "IGM"
    if os.path.exists(igm_dir):
        for item in os.listdir(igm_dir):
            s = igm_dir / item
            d = src_dir / "igm" / item
            copy_item(s, d)
    
    # Move files from other directories
    for old_dir, new_subdir in [
        ("catalog", "catalog"), 
        ("halo", "halo"), 
        ("lensing", "lensing"), 
        ("rbstat", "stats"), 
        ("utils", "utils")
    ]:
        print(f"Moving {old_dir} files...")
        source_dir = root_dir / old_dir
        if os.path.exists(source_dir):
            for item in os.listdir(source_dir):
                s = source_dir / item
                d = src_dir / new_subdir / item
                copy_item(s, d)
    
    # Copy files from tests directory, being careful not to copy to self
    print("Moving test files...")
    test_dir = root_dir / "tests"
    if os.path.exists(test_dir):
        for item in os.listdir(test_dir):
            s = test_dir / item
            
            # Skip copying files to themselves
            if s.resolve() == (test_dir / item).resolve():
                if item == "__init__.py" or not os.path.isdir(s):
                    # Special case for __init__.py and other files: 
                    # don't try to move them within the same directory
                    continue
            
            # Try to determine which module this test belongs to
            if item.startswith("test_") and "_" in item and os.path.isfile(s):
                module_name = item.split("_")[1].split(".")[0]  # Extract module name from test file
                if module_name in subdirs:
                    d = root_dir / "tests" / f"test_{module_name}" / item
                    copy_item(s, d)
                    continue
            
            # For other files/dirs, keep them at top level of tests
            d = root_dir / "tests" / item
            # Only copy if source and destination are different
            if s.resolve() != d.resolve():
                copy_item(s, d)
    
    # Copy docs
    print("Moving doc files...")
    doc_dir = root_dir / "docs"
    if os.path.exists(doc_dir):
        for item in os.listdir(doc_dir):
            s = doc_dir / item
            if item == "_build":
                # Skip _build directory to avoid copying built docs
                continue
            d = root_dir / "docs" / item
            # Only copy if source and destination are different
            if s.resolve() != d.resolve():
                copy_item(s, d)
    
    # Create init files
    print("Creating __init__.py files...")
    for path in [src_dir] + [src_dir / subdir for subdir in subdirs]:
        init_file = path / "__init__.py"
        if not os.path.exists(init_file):
            with open(init_file, "w") as f:
                if path == src_dir:
                    f.write('"""rbcodes: A package for analyzing astronomical spectroscopic data."""\n\n')
                    f.write('__version__ = "0.1.0"\n')
                else:
                    f.write('"""rbcodes {} module."""\n'.format(path.name))
            print(f"Created {init_file}")
    
    # Create setup.py
    if not os.path.exists(root_dir / "setup.py"):
        setup_py = """
from setuptools import setup, find_packages

setup(
    name="rbcodes",
    version="0.1.0",
    description="A package for analyzing astronomical spectroscopic data",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/rbcodes",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        # Add your dependencies here
    ],
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
        """
        with open(root_dir / "setup.py", "w") as f:
            f.write(setup_py.strip())
        print("Created setup.py")
    
    # Create pyproject.toml
    if not os.path.exists(root_dir / "pyproject.toml"):
        pyproject_toml = """
[build-system]
requires = ["setuptools>=42", "wheel"]
build-backend = "setuptools.build_meta"

[tool.black]
line-length = 88
target-version = ['py37', 'py38', 'py39']
include = '\.pyi?$'

[tool.isort]
profile = "black"
multi_line_output = 3
        """
        with open(root_dir / "pyproject.toml", "w") as f:
            f.write(pyproject_toml.strip())
        print("Created pyproject.toml")
    
    # Create scripts for running tools
    for script_name, module_path in [
        ("run_abstools.py", "rbcodes.gui.abstools.AbsToolsLauncher"),
        ("run_multispecviewer.py", "rbcodes.gui.multispecviewer.multispec"),
        ("run_spectral_gui.py", "rbcodes.gui.spectral_gui.run_spectral_gui"),
    ]:
        script_path = root_dir / "scripts" / script_name
        if not os.path.exists(script_path):
            script_content = f"""#!/usr/bin/env python3
import sys
from {module_path} import main

if __name__ == "__main__":
    sys.exit(main())
"""
            with open(script_path, "w") as f:
                f.write(script_content)
            os.chmod(script_path, 0o755)  # Make executable
            print(f"Created {script_name}")

    print("\nRestructuring complete!")
    print("Note: This script copies files but doesn't remove the originals.")
    print("After verifying everything works correctly, you may want to remove the original directories.")
    print("Make sure to update imports in your code to use the new module structure.")

if __name__ == "__main__":
    main()