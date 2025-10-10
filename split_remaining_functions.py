#!/usr/bin/env python3

"""
Script to complete the splitting of functions.R into multiple files
This script extracts all remaining functions and splits them into files with max 5 functions each
"""

import re
import os
import glob
from pathlib import Path
import math

def read_file(filepath):
    """Read the entire content of a file"""
    with open(filepath, 'r', encoding='utf-8') as f:
        return f.read()

def write_file(content, filepath):
    """Write content to a file"""
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

def extract_functions(content):
    """Extract all R function definitions from the content"""
    # Pattern to match R function definitions: name <- function(params) {
    function_pattern = r'([a-zA-Z_][a-zA-Z0-9_.]*)\s*<-\s*function\s*\([^)]*\)\s*\{'
    
    # Find all function matches
    function_matches = list(re.finditer(function_pattern, content))
    
    if not function_matches:
        print("No functions found.")
        return {}
    
    # Extract function names
    function_names = [match.group(1) for match in function_matches]
    
    # Find start and end positions for each function
    functions_dict = {}
    
    for i, match in enumerate(function_matches):
        start_pos = match.start()
        func_name = function_names[i]
        
        # Find the end of the function by counting braces
        brace_count = 0
        pos = start_pos
        in_function_def = False
        end_pos = -1
        
        while pos < len(content):
            char = content[pos]
            
            if char == '{':
                brace_count += 1
                in_function_def = True
            elif char == '}':
                brace_count -= 1
                if in_function_def and brace_count == 0:
                    end_pos = pos + 1  # Include the closing brace
                    break
            pos += 1
        
        # Extract the complete function text
        if end_pos > 0:
            func_text = content[start_pos:end_pos]
            functions_dict[func_name] = func_text
    
    return functions_dict

def get_file_theme(func_names, functions_dict):
    """Determine the theme for a group of functions based on their content"""
    # Combine all function text for analysis
    func_text = ' '.join([functions_dict.get(name, '') for name in func_names])
    func_text = func_text.lower()
    
    if re.search(r'plot|ggplot|featureplot|dimplot', func_text):
        return "plot_functions"
    elif re.search(r'diffex|findmarkers|enrichment', func_text):
        return "diffex_functions"
    elif re.search(r'numbat|clone|heatmap', func_text):
        return "numbat_functions"
    elif re.search(r'cluster|seurat|seu', func_text):
        return "cluster_functions"
    elif re.search(r'score|pca|qc|filter', func_text):
        return "scoring_functions"
    elif re.search(r'read|write|save|load', func_text):
        return "io_functions"
    else:
        return "utility_functions"

def get_existing_files():
    """Get list of existing R files in the R/ directory"""
    r_dir = Path("R")
    if not r_dir.exists():
        r_dir.mkdir(exist_ok=True)
    
    return list(r_dir.glob("*.R"))

def get_next_file_number(theme, existing_files):
    """Get the next available file number for a given theme"""
    theme_files = [f for f in existing_files if f.name.startswith(theme)]
    
    if not theme_files:
        return 1
    
    # Extract numbers from theme files
    numbers = []
    for f in theme_files:
        match = re.search(r'(\d+)\.R$', f.name)
        if match:
            numbers.append(int(match.group(1)))
    
    return max(numbers) + 1 if numbers else 1

def split_functions_into_groups(functions_dict, functions_per_file=5):
    """Split functions into groups of specified size"""
    func_names = list(functions_dict.keys())
    groups = []
    
    for i in range(0, len(func_names), functions_per_file):
        groups.append(func_names[i:i + functions_per_file])
    
    return groups

def create_file_header(theme, file_number):
    """Create a header for the function file"""
    theme_title = theme.replace('_', ' ').title()
    return f"# {theme_title} ({file_number})\n\n"

def main():
    """Main function to split the functions"""
    # Check if functions.R exists
    functions_file = Path("functions.R")
    if not functions_file.exists():
        print("Error: functions.R file not found!")
        return
    
    # Read the original functions.R file
    print("Reading functions.R...")
    functions_content = read_file(functions_file)
    
    # Extract all functions
    print("Extracting functions...")
    functions_dict = extract_functions(functions_content)
    
    if not functions_dict:
        print("No functions found to split.")
        return
    
    print(f"Found {len(functions_dict)} functions:")
    for name in list(functions_dict.keys())[:10]:  # Show first 10
        print(f"  - {name}")
    if len(functions_dict) > 10:
        print(f"  ... and {len(functions_dict) - 10} more")
    
    # Group functions into files of 5 each
    function_groups = split_functions_into_groups(functions_dict, 5)
    
    # Get existing files to avoid conflicts
    existing_files = get_existing_files()
    existing_count = len(existing_files)
    
    print(f"\nSplitting into {len(function_groups)} files...")
    
    # Create files for each function group
    created_files = 0
    for i, group_names in enumerate(function_groups):
        theme = get_file_theme(group_names, functions_dict)
        
        # Get next available file number for this theme
        file_number = get_next_file_number(theme, existing_files)
        filename = f"R/{theme}_{file_number}.R"
        
        # If file already exists, use a different numbering scheme
        if Path(filename).exists():
            filename = f"R/{theme}_{existing_count + i + 1}.R"
        
        # Create file content
        file_content = create_file_header(theme, file_number)
        
        # Add each function to the file
        for func_name in group_names:
            if func_name in functions_dict:
                file_content += functions_dict[func_name] + "\n\n"
        
        # Write the file
        write_file(file_content, filename)
        created_files += 1
        
        print(f"Created: {filename} with functions: {', '.join(group_names)}")
        
        # Update existing files list to avoid conflicts
        existing_files.append(Path(filename))
    
    print(f"\nFunction splitting completed!")
    print(f"Total files created: {created_files}")
    print("Use source('R/load_all_functions.R') to load all functions.")

if __name__ == "__main__":
    main()