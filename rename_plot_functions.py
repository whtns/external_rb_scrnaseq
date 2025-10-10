#!/usr/bin/env python3

"""
Script to rename plot_functions files in sequential order
This script will rename all plot_functions_*.R files to have sequential numbering (1, 2, 3, ...)
"""

import os
import re
from pathlib import Path
import shutil

def get_plot_function_files():
    """Get all plot_functions files sorted by their current number"""
    r_dir = Path("R")
    if not r_dir.exists():
        print("R/ directory not found!")
        return []
    
    # Find all plot_functions files
    plot_files = list(r_dir.glob("plot_functions_*.R"))
    
    # Extract numbers and sort
    def extract_number(filename):
        match = re.search(r'plot_functions_(\d+)\.R$', filename.name)
        return int(match.group(1)) if match else 0
    
    plot_files.sort(key=extract_number)
    return plot_files

def create_backup():
    """Create a backup directory for safety"""
    backup_dir = Path("R_backup_rename")
    backup_dir.mkdir(exist_ok=True)
    return backup_dir

def update_load_all_functions(old_to_new_mapping):
    """Update the load_all_functions.R file with new filenames"""
    load_script = Path("R/load_all_functions.R")
    if not load_script.exists():
        return
    
    try:
        content = load_script.read_text()
        
        # Replace old filenames with new ones in the load script
        for old_name, new_name in old_to_new_mapping.items():
            old_pattern = f'source("R/{old_name}")'
            new_pattern = f'source("R/{new_name}")'
            content = content.replace(old_pattern, new_pattern)
        
        load_script.write_text(content)
        print(f"Updated {load_script.name}")
        
    except Exception as e:
        print(f"Error updating load_all_functions.R: {e}")

def main():
    """Main function to rename plot_functions files sequentially"""
    print("=== Renaming plot_functions files in sequential order ===\n")
    
    # Get all plot_functions files
    plot_files = get_plot_function_files()
    
    if not plot_files:
        print("No plot_functions files found!")
        return
    
    print(f"Found {len(plot_files)} plot_functions files")
    print("Current files:")
    for i, file in enumerate(plot_files, 1):
        print(f"  {i:2d}. {file.name}")
    
    # Create backup
    backup_dir = create_backup()
    print(f"\nCreating backups in {backup_dir}/")
    
    # Create mapping from old names to new names
    old_to_new_mapping = {}
    temp_names = []
    
    # First, rename all files to temporary names to avoid conflicts
    print("\nStep 1: Renaming to temporary names...")
    for i, old_file in enumerate(plot_files, 1):
        # Create backup
        backup_path = backup_dir / old_file.name
        shutil.copy2(old_file, backup_path)
        
        # Create temporary name
        temp_name = f"plot_functions_temp_{i}.R"
        temp_path = old_file.parent / temp_name
        
        # Rename to temporary name
        old_file.rename(temp_path)
        temp_names.append((temp_path, i))
        
        print(f"  {old_file.name} -> {temp_name}")
    
    # Second, rename temporary files to final sequential names
    print("\nStep 2: Renaming to final sequential names...")
    for temp_path, new_number in temp_names:
        new_name = f"plot_functions_{new_number}.R"
        new_path = temp_path.parent / new_name
        
        # Extract original name for mapping
        original_name = backup_dir.glob(f"*temp_{new_number}.R")
        # Find the original name from backup directory
        for backup_file in backup_dir.iterdir():
            if backup_file.name.startswith("plot_functions_") and backup_file.name.endswith(".R"):
                # Read both files and compare content to find the match
                if backup_file.stat().st_size == temp_path.stat().st_size:
                    old_to_new_mapping[backup_file.name] = new_name
                    break
        
        # Rename to final name
        temp_path.rename(new_path)
        print(f"  {temp_path.name} -> {new_name}")
    
    # Create proper mapping by comparing file sizes and positions
    print("\nCreating filename mapping...")
    backup_files = sorted([f for f in backup_dir.iterdir() if f.name.startswith("plot_functions_")], 
                         key=lambda x: int(re.search(r'plot_functions_(\d+)\.R$', x.name).group(1)))
    
    old_to_new_mapping = {}
    for i, backup_file in enumerate(backup_files, 1):
        new_name = f"plot_functions_{i}.R"
        old_to_new_mapping[backup_file.name] = new_name
        print(f"  {backup_file.name} -> {new_name}")
    
    # Update load_all_functions.R if it exists
    print(f"\nUpdating references...")
    update_load_all_functions(old_to_new_mapping)
    
    # Final verification
    print(f"\nVerification - Final files:")
    final_files = sorted(Path("R").glob("plot_functions_*.R"), 
                        key=lambda x: int(re.search(r'plot_functions_(\d+)\.R$', x.name).group(1)))
    
    for i, file in enumerate(final_files, 1):
        expected_name = f"plot_functions_{i}.R"
        if file.name == expected_name:
            print(f"  ✓ {file.name}")
        else:
            print(f"  ✗ {file.name} (expected: {expected_name})")
    
    print(f"\n=== Summary ===")
    print(f"Total files renamed: {len(plot_files)}")
    print(f"New sequential range: plot_functions_1.R to plot_functions_{len(plot_files)}.R")
    print(f"Backups saved in: {backup_dir}/")
    print("Renaming completed successfully!")

if __name__ == "__main__":
    main()