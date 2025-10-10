#!/usr/bin/env python3

"""
Script to add roxygen2 documentation to all R functions
This script analyzes R functions and adds appropriate roxygen2 documentation
based on function names, parameters, and content analysis.
"""

import re
import os
from pathlib import Path
from typing import Dict, List, Tuple
import ast

def read_file(filepath):
    """Read the entire content of a file"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            return f.read()
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return ""

def write_file(content, filepath):
    """Write content to a file"""
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

def extract_function_info(content):
    """Extract function information including parameters and body"""
    # Pattern to match R function definitions
    function_pattern = r'([a-zA-Z_][a-zA-Z0-9_.]*)\s*<-\s*function\s*\(([^)]*)\)\s*\{'
    
    functions = []
    function_matches = list(re.finditer(function_pattern, content))
    
    for i, match in enumerate(function_matches):
        start_pos = match.start()
        func_name = match.group(1)
        params_str = match.group(2).strip()
        
        # Parse parameters
        params = parse_function_parameters(params_str)
        
        # Find the end of the function
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
                    end_pos = pos + 1
                    break
            pos += 1
        
        if end_pos > 0:
            func_text = content[start_pos:end_pos]
            
            # Find the line where the function starts
            lines_before = content[:start_pos].count('\n')
            
            functions.append({
                'name': func_name,
                'params': params,
                'text': func_text,
                'start_pos': start_pos,
                'end_pos': end_pos,
                'start_line': lines_before,
                'full_params_str': params_str
            })
    
    return functions

def parse_function_parameters(params_str):
    """Parse R function parameters string"""
    if not params_str.strip():
        return []
    
    params = []
    # Split by comma, but be careful of nested structures
    param_parts = []
    paren_count = 0
    current_part = ""
    
    for char in params_str:
        if char == ',' and paren_count == 0:
            param_parts.append(current_part.strip())
            current_part = ""
        else:
            if char in '([{':
                paren_count += 1
            elif char in ')]}':
                paren_count -= 1
            current_part += char
    
    if current_part.strip():
        param_parts.append(current_part.strip())
    
    for part in param_parts:
        if '=' in part:
            name, default = part.split('=', 1)
            params.append({
                'name': name.strip(),
                'default': default.strip(),
                'has_default': True
            })
        else:
            params.append({
                'name': part.strip(),
                'default': None,
                'has_default': False
            })
    
    return params

def generate_function_description(func_name, func_text):
    """Generate a description for the function based on its name and content"""
    name_lower = func_name.lower()
    
    # Common patterns in function names
    if 'plot' in name_lower:
        if 'numbat' in name_lower:
            return f"Create a numbat-related plot visualization"
        elif 'volcano' in name_lower:
            return f"Generate a volcano plot for differential expression data"
        elif 'heatmap' in name_lower:
            return f"Create a heatmap visualization"
        elif 'feature' in name_lower:
            return f"Generate a feature plot"
        else:
            return f"Create a plot visualization"
    
    elif 'filter' in name_lower:
        return f"Filter data based on specified criteria"
    
    elif 'score' in name_lower or 'scoring' in name_lower:
        return f"Calculate scores for the given data"
    
    elif 'diffex' in name_lower or 'differential' in name_lower:
        return f"Perform differential expression analysis"
    
    elif 'enrich' in name_lower:
        return f"Perform enrichment analysis"
    
    elif 'numbat' in name_lower:
        return f"Perform numbat-related analysis"
    
    elif 'cluster' in name_lower:
        return f"Perform clustering analysis"
    
    elif 'load' in name_lower or 'read' in name_lower:
        return f"Load or read data from file"
    
    elif 'save' in name_lower or 'write' in name_lower:
        return f"Save or write data to file"
    
    elif 'convert' in name_lower:
        return f"Convert data from one format to another"
    
    elif 'extract' in name_lower or 'pull' in name_lower:
        return f"Extract or pull specific data elements"
    
    elif 'merge' in name_lower or 'join' in name_lower:
        return f"Merge or join datasets"
    
    elif 'calculate' in name_lower or 'compute' in name_lower:
        return f"Calculate or compute values"
    
    elif 'annotate' in name_lower:
        return f"Add annotations to data"
    
    elif 'normalize' in name_lower:
        return f"Normalize data values"
    
    else:
        # Generic description
        return f"Perform {func_name.replace('_', ' ').lower()} operation"

def infer_parameter_description(param_name, param_default, func_text):
    """Infer parameter description based on name and usage"""
    name_lower = param_name.lower()
    
    # Common parameter patterns
    if param_name == '...':
        return "Additional arguments passed to other functions"
    
    elif name_lower in ['data', 'df', 'dataset']:
        return "Input data frame or dataset"
    
    elif name_lower in ['seu', 'seurat', 'myseu']:
        return "Seurat object"
    
    elif name_lower in ['nb', 'numbat', 'mynb']:
        return "Numbat object"
    
    elif name_lower in ['plot', 'p', 'myplot']:
        return "Plot object (ggplot2)"
    
    elif 'file' in name_lower or 'path' in name_lower:
        return "File path"
    
    elif 'title' in name_lower:
        return "Plot title"
    
    elif 'ident' in name_lower:
        return "Cell identities or groups"
    
    elif 'cluster' in name_lower:
        return "Cluster information"
    
    elif 'gene' in name_lower:
        return "Gene names or identifiers"
    
    elif 'cell' in name_lower:
        return "Cell identifiers or information"
    
    elif 'threshold' in name_lower or 'cutoff' in name_lower:
        return "Threshold value for filtering"
    
    elif 'color' in name_lower or 'col' in name_lower:
        return "Color specification"
    
    elif 'size' in name_lower:
        return "Size parameter"
    
    elif 'alpha' in name_lower:
        return "Transparency level (0-1)"
    
    elif name_lower in ['method', 'type']:
        return "Method or type specification"
    
    elif name_lower in ['verbose', 'quiet']:
        return "Whether to print progress messages"
    
    elif param_default == 'TRUE' or param_default == 'FALSE':
        return f"Logical flag (default: {param_default})"
    
    elif param_default and param_default.startswith('"'):
        return f"Character string (default: {param_default})"
    
    else:
        return f"Parameter for {param_name.replace('_', ' ')}"

def infer_return_value(func_name, func_text):
    """Infer what the function returns based on its content and name"""
    name_lower = func_name.lower()
    
    # Look for return statements
    if 'return(' in func_text:
        # Try to extract what's being returned
        return_match = re.search(r'return\s*\(\s*([^)]+)\s*\)', func_text)
        if return_match:
            returned_obj = return_match.group(1).strip()
            if 'plot' in returned_obj.lower() or 'ggplot' in func_text:
                return "ggplot2 plot object"
            elif 'seu' in returned_obj.lower():
                return "Modified Seurat object"
            elif 'data' in returned_obj.lower() or 'df' in returned_obj.lower():
                return "Data frame"
            elif 'list' in returned_obj.lower():
                return "List object"
    
    # Infer from function name
    if 'plot' in name_lower:
        return "ggplot2 plot object"
    elif 'filter' in name_lower:
        return "Filtered data"
    elif 'score' in name_lower:
        return "Numeric scores or scored data"
    elif 'diffex' in name_lower:
        return "Differential expression results"
    elif 'enrich' in name_lower:
        return "Enrichment analysis results"
    elif 'load' in name_lower or 'read' in name_lower:
        return "Loaded data object"
    elif 'convert' in name_lower:
        return "Converted data object"
    elif 'extract' in name_lower or 'pull' in name_lower:
        return "Extracted data elements"
    elif 'calculate' in name_lower or 'compute' in name_lower:
        return "Calculated values"
    else:
        return "Function result"

def generate_roxygen_docs(func_info):
    """Generate roxygen2 documentation for a function"""
    func_name = func_info['name']
    params = func_info['params']
    func_text = func_info['text']
    
    # Generate description
    description = generate_function_description(func_name, func_text)
    
    # Generate parameter documentation
    param_docs = []
    for param in params:
        param_desc = infer_parameter_description(param['name'], param.get('default'), func_text)
        param_docs.append(f"#' @param {param['name']} {param_desc}")
    
    # Generate return documentation
    return_desc = infer_return_value(func_name, func_text)
    
    # Construct full documentation
    docs = [
        f"#' {description}",
        "#'"
    ]
    
    # Add parameter docs
    docs.extend(param_docs)
    
    # Add return documentation
    docs.append(f"#' @return {return_desc}")
    docs.append("#' @export")
    
    return '\n'.join(docs) + '\n'

def add_roxygen_to_file(filepath):
    """Add roxygen documentation to all functions in a file"""
    print(f"Processing {filepath.name}...")
    
    content = read_file(filepath)
    if not content:
        return False
    
    # Extract all functions
    functions = extract_function_info(content)
    if not functions:
        print(f"  No functions found in {filepath.name}")
        return False
    
    print(f"  Found {len(functions)} functions")
    
    # Add documentation in reverse order to preserve positions
    functions.sort(key=lambda x: x['start_pos'], reverse=True)
    
    modified_content = content
    
    for func_info in functions:
        # Check if function already has roxygen documentation
        lines = modified_content[:func_info['start_pos']].split('\n')
        prev_lines = lines[-5:] if len(lines) >= 5 else lines
        
        has_roxygen = any(line.strip().startswith("#'") for line in prev_lines)
        
        if has_roxygen:
            print(f"    Skipping {func_info['name']} - already has roxygen docs")
            continue
        
        # Generate roxygen documentation
        roxygen_docs = generate_roxygen_docs(func_info)
        
        # Insert documentation before function
        before = modified_content[:func_info['start_pos']]
        after = modified_content[func_info['start_pos']:]
        
        # Make sure there's proper spacing
        if before and not before.endswith('\n'):
            before += '\n'
        
        modified_content = before + roxygen_docs + after
        
        print(f"    Added docs for {func_info['name']}")
    
    # Write the modified content
    write_file(modified_content, filepath)
    return True

def main():
    """Main function to add roxygen documentation to all R files"""
    print("=== Adding Roxygen2 Documentation to R Functions ===\n")
    
    r_dir = Path("R")
    if not r_dir.exists():
        print("R/ directory not found!")
        return
    
    # Get all R files except the original functions.R and load script
    r_files = [f for f in r_dir.glob("*.R") 
               if f.name not in ['functions.R', 'load_all_functions.R']]
    
    if not r_files:
        print("No R function files found!")
        return
    
    print(f"Found {len(r_files)} R files to process\n")
    
    # Create backup directory
    backup_dir = Path("R_backup_roxygen")
    backup_dir.mkdir(exist_ok=True)
    print(f"Creating backups in {backup_dir}/\n")
    
    processed_files = 0
    total_functions = 0
    
    for filepath in sorted(r_files):
        # Create backup
        backup_path = backup_dir / filepath.name
        content = read_file(filepath)
        write_file(content, backup_path)
        
        # Add roxygen documentation
        if add_roxygen_to_file(filepath):
            processed_files += 1
            # Count functions for summary
            functions = extract_function_info(read_file(filepath))
            total_functions += len(functions)
        
        print()  # Empty line for readability
    
    print("=== Summary ===")
    print(f"Files processed: {processed_files}")
    print(f"Total functions documented: {total_functions}")
    print(f"Backups saved in: {backup_dir}/")
    print("\nRoxygen documentation added successfully!")
    print("\nTo generate documentation, run:")
    print("  devtools::document()")
    print("or")
    print("  roxygen2::roxygenise()")

if __name__ == "__main__":
    main()