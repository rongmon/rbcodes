#!/usr/bin/env python3
"""
Docstring Extractor

This script extracts docstrings from Python files in a specified directory
and generates a README.md file with the documentation.

chmod +x docstring_extractor.py

Usage:
    python docstring_extractor.py [directory] [--output README.md]
    # Basic usage (current directory, outputs to README.md)
    python docstring_extractor.py

    # Specify a directory
    python docstring_extractor.py path/to/your/project

    # Specify an output file
    python docstring_extractor.py --output docs/API.md

    # Both options
    python docstring_extractor.py path/to/your/project --output custom_readme.md
"""

import os
import sys
import ast
import re
import argparse
from typing import Dict, List, Tuple, Any, Optional


class DocstringVisitor(ast.NodeVisitor):
    """Visit AST nodes and extract docstrings from modules, classes, and functions."""
    
    def __init__(self):
        self.docstrings = {
            'module': {},
            'class': {},
            'function': {},
        }
        self.current_module = None
    
    def visit_Module(self, node: ast.Module):
        """Extract module docstring."""
        if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Str):
            self.docstrings['module'][self.current_module] = node.body[0].value.s.strip()
        self.generic_visit(node)
    
    def visit_ClassDef(self, node: ast.ClassDef):
        """Extract class docstring."""
        if (node.body and isinstance(node.body[0], ast.Expr) and 
                isinstance(node.body[0].value, ast.Str)):
            class_key = f"{self.current_module}.{node.name}"
            self.docstrings['class'][class_key] = {
                'name': node.name,
                'docstring': node.body[0].value.s.strip()
            }
        self.generic_visit(node)
    
    def visit_FunctionDef(self, node: ast.FunctionDef):
        """Extract function docstring."""
        if (node.body and isinstance(node.body[0], ast.Expr) and 
                isinstance(node.body[0].value, ast.Str)):
            function_key = f"{self.current_module}.{node.name}"
            self.docstrings['function'][function_key] = {
                'name': node.name,
                'docstring': node.body[0].value.s.strip()
            }
        self.generic_visit(node)


def extract_docstrings(directory: str) -> Dict:
    """
    Extract docstrings from all Python files in the given directory.
    
    Args:
        directory: Path to the directory containing Python files
        
    Returns:
        A dictionary with extracted docstrings
    """
    visitor = DocstringVisitor()
    
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                file_path = os.path.join(root, file)
                rel_path = os.path.relpath(file_path, directory)
                module_name = os.path.splitext(rel_path)[0].replace(os.path.sep, '.')
                
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        code = f.read()
                    
                    visitor.current_module = module_name
                    tree = ast.parse(code)
                    visitor.visit(tree)
                except Exception as e:
                    print(f"Error parsing {file_path}: {e}")
    
    return visitor.docstrings


def generate_readme(docstrings: Dict, output_file: str = 'README.md') -> None:
    """
    Generate a README.md file from the extracted docstrings.
    
    Args:
        docstrings: Dictionary with extracted docstrings
        output_file: Path to the output README.md file
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        # Write main title
        f.write("# Project Documentation\n\n")
        f.write("*Auto-generated documentation from docstrings*\n\n")
        
        # Write module docstrings
        if docstrings['module']:
            f.write("## Modules\n\n")
            for module_name, docstring in docstrings['module'].items():
                f.write(f"### {module_name}\n\n")
                f.write(f"{docstring}\n\n")
        
        # Write class docstrings
        if docstrings['class']:
            f.write("## Classes\n\n")
            for class_key, info in docstrings['class'].items():
                module_name, class_name = class_key.rsplit('.', 1)
                f.write(f"### {class_name} (`{module_name}`)\n\n")
                f.write(f"{info['docstring']}\n\n")
        
        # Write function docstrings
        if docstrings['function']:
            f.write("## Functions\n\n")
            for func_key, info in docstrings['function'].items():
                module_name, func_name = func_key.rsplit('.', 1)
                f.write(f"### {func_name}() (`{module_name}`)\n\n")
                f.write(f"{info['docstring']}\n\n")


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description='Extract docstrings and generate README.md')
    parser.add_argument('directory', nargs='?', default='.', 
                        help='Directory containing Python files (default: current directory)')
    parser.add_argument('--output', '-o', default='README.md',
                        help='Output file (default: README.md)')
    
    args = parser.parse_args()
    
    print(f"Extracting docstrings from {args.directory}...")
    docstrings = extract_docstrings(args.directory)
    
    print(f"Generating {args.output}...")
    generate_readme(docstrings, args.output)
    
    print("Done!")


if __name__ == "__main__":
    main()
