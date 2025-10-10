#!/usr/bin/env python3

"""
Script to identify and remove duplicated function definitions across R files
This script analyzes all R files in the R/ directory, identifies duplicate functions,
and removes them while preserving one copy of each unique function.
"""

import re
import os
from pathlib import Path
from collections import defaultdict
import hashlib

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

def extract_functions_from_content(content, filepath):
    """Extract all R function definitions from content with their positions"""
    # Pattern to match R function definitions: name <- function(params) {
    function_pattern = r'([a-zA-Z_][a-zA-Z0-9_.]*)\s*<-\s*function\s*\([^)]*\)\s*\{'
    
    functions = []
    function_matches = list(re.finditer(function_pattern, content))
    
    for i, match in enumerate(function_matches):
        start_pos = match.start()
        func_name = match.group(1)
        
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
        
        if end_pos > 0:
            func_text = content[start_pos:end_pos]
            # Normalize whitespace for comparison
            normalized_text = re.sub(r'\s+', ' ', func_text.strip())
            
            functions.append({
                'name': func_name,
                'text': func_text,
                'normalized': normalized_text,
                'start': start_pos,
                'end': end_pos,
                'file': filepath,
                'hash': hashlib.md5(normalized_text.encode()).hexdigest()
            })
    
    return functions

def analyze_all_functions():
    """Analyze all R files and extract function information"""
    r_dir = Path("R")
    if not r_dir.exists():
        print("R/ directory not found!")
        return {}, {}
    
    all_functions = []
    file_contents = {}
    
    # Get all R files except the original functions.R and load script
    r_files = [f for f in r_dir.glob("*.R") 
               if f.name not in ['functions.R', 'load_all_functions.R']]
    
    print(f"Analyzing {len(r_files)} R files...")
    
    for filepath in r_files:
        content = read_file(filepath)
        if content:
            file_contents[filepath] = content
            functions = extract_functions_from_content(content, filepath)
            all_functions.extend(functions)
            if functions:
                print(f"  {filepath.name}: {len(functions)} functions")
    
    return all_functions, file_contents

def find_duplicates(all_functions):
    """Find duplicate functions based on normalized content"""
    # Group by function name and content hash
    name_groups = defaultdict(list)
    hash_groups = defaultdict(list)
    
    for func in all_functions:
        name_groups[func['name']].append(func)
        hash_groups[func['hash']].append(func)
    
    duplicates_by_name = {name: funcs for name, funcs in name_groups.items() if len(funcs) > 1}
    duplicates_by_content = {hash_val: funcs for hash_val, funcs in hash_groups.items() if len(funcs) > 1}
    
    return duplicates_by_name, duplicates_by_content

def remove_duplicates_from_file(filepath, functions_to_remove, content):
    """Remove specified functions from a file's content"""
    # Sort by position in reverse order to avoid offset issues
    functions_to_remove.sort(key=lambda x: x['start'], reverse=True)
    
    modified_content = content
    removed_count = 0
    
    for func in functions_to_remove:
        # Remove the function from the content
        before = modified_content[:func['start']]
        after = modified_content[func['end']:]
        
        # Clean up extra newlines
        before = before.rstrip()
        after = after.lstrip('\n')
        if after and not after.startswith('\n'):
            after = '\n' + after
        
        modified_content = before + after
        removed_count += 1
    
    return modified_content, removed_count

def create_backup_directory():
    """Create a backup directory for original files"""
    backup_dir = Path("R_backup")
    backup_dir.mkdir(exist_ok=True)
    return backup_dir

def main():
    """Main function to remove duplicate functions"""
    print("=== R Function Deduplication Script ===\n")
    
    # Analyze all functions
    all_functions, file_contents = analyze_all_functions()
    
    if not all_functions:
        print("No functions found!")
        return
    
    print(f"\nTotal functions found: {len(all_functions)}")
    
    # Find duplicates
    duplicates_by_name, duplicates_by_content = find_duplicates(all_functions)
    
    print(f"Functions with same name: {len(duplicates_by_name)}")
    print(f"Functions with identical content: {len(duplicates_by_content)}")
    
    if not duplicates_by_name and not duplicates_by_content:
        print("No duplicates found!")
        return
    
    # Create backup directory
    backup_dir = create_backup_directory()
    print(f"\nCreating backups in {backup_dir}/...")
    
    # Strategy: For each set of duplicates, keep the first occurrence and remove the rest
    files_to_modify = {}  # filepath -> list of functions to remove
    
    # Handle duplicates by content (exact duplicates)
    print("\nAnalyzing content duplicates...")
    for hash_val, duplicate_funcs in duplicates_by_content.items():
        if len(duplicate_funcs) <= 1:
            continue
            
        # Sort by file creation time or name to be consistent
        duplicate_funcs.sort(key=lambda x: (x['file'].name, x['start']))
        keep_func = duplicate_funcs[0]
        remove_funcs = duplicate_funcs[1:]
        
        print(f"  Function '{keep_func['name']}' found {len(duplicate_funcs)} times")
        print(f"    Keeping: {keep_func['file'].name}")
        for func in remove_funcs:
            print(f"    Removing from: {func['file'].name}")
            if func['file'] not in files_to_modify:
                files_to_modify[func['file']] = []
            files_to_modify[func['file']].append(func)
    
    # Handle name duplicates with different content (potential issues)
    print("\nAnalyzing name duplicates with different content...")
    for name, duplicate_funcs in duplicates_by_name.items():
        if len(duplicate_funcs) <= 1:
            continue
            
        # Check if they have the same content hash
        hashes = set(func['hash'] for func in duplicate_funcs)
        if len(hashes) == 1:
            continue  # Already handled in content duplicates
        
        print(f"  WARNING: Function '{name}' has different implementations:")
        for func in duplicate_funcs:
            print(f"    In {func['file'].name} (hash: {func['hash'][:8]}...)")
        print(f"    Keeping first occurrence in {duplicate_funcs[0]['file'].name}")
        
        # Remove all but the first
        for func in duplicate_funcs[1:]:
            if func['file'] not in files_to_modify:
                files_to_modify[func['file']] = []
            files_to_modify[func['file']].append(func)
    
    if not files_to_modify:
        print("No files need modification!")
        return
    
    # Backup and modify files
    total_removed = 0
    empty_files = []
    
    for filepath, functions_to_remove in files_to_modify.items():
        # Create backup
        backup_path = backup_dir / filepath.name
        write_file(file_contents[filepath], backup_path)
        
        # Remove duplicates
        modified_content, removed_count = remove_duplicates_from_file(
            filepath, functions_to_remove, file_contents[filepath]
        )
        
        # Check if file is now empty (except for header comments)
        content_without_comments = re.sub(r'^\s*#.*$', '', modified_content, flags=re.MULTILINE).strip()
        
        if not content_without_comments:
            empty_files.append(filepath)
            print(f"  {filepath.name}: Removed {removed_count} functions - FILE NOW EMPTY")
        else:
            # Write modified content
            write_file(modified_content, filepath)
            print(f"  {filepath.name}: Removed {removed_count} functions")
        
        total_removed += removed_count
    
    # Handle empty files
    if empty_files:
        print(f"\nFound {len(empty_files)} empty files:")
        for filepath in empty_files:
            print(f"  - {filepath.name}")
        
        response = input("\nDelete empty files? (y/N): ").strip().lower()
        if response == 'y':
            for filepath in empty_files:
                filepath.unlink()
                print(f"  Deleted: {filepath.name}")
    
    # Summary
    print(f"\n=== Summary ===")
    print(f"Total duplicate functions removed: {total_removed}")
    print(f"Files modified: {len(files_to_modify)}")
    print(f"Backups created in: {backup_dir}/")
    
    # Update load_all_functions.R if it exists
    load_script = Path("R/load_all_functions.R")
    if load_script.exists():
        print("\nUpdating load_all_functions.R...")
        remaining_files = [f for f in Path("R").glob("*.R") 
                          if f.name not in ['functions.R', 'load_all_functions.R'] and f.exists()]
        
        load_content = "# Load all function files\n\n"
        for filepath in sorted(remaining_files):
            load_content += f'source("R/{filepath.name}")\n'
        
        write_file(load_content, load_script)
        print(f"Updated {load_script.name} with {len(remaining_files)} files")
    
    print("\nDeduplication complete!")

if __name__ == "__main__":
    main()