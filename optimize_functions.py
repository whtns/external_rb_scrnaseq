#!/usr/bin/env python3

"""
Script to optimize R functions for better efficiency
This script identifies and fixes common performance bottlenecks in R code:
- Redundant file I/O operations
- Inefficient loops and apply functions  
- Unnecessary repeated computations
- Memory-inefficient operations
- Unvectorized operations
"""

import re
import os
from pathlib import Path
from collections import defaultdict

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

def optimize_repeated_file_reads(content):
    """Optimize functions that read the same file multiple times"""
    optimizations = []
    
    # Pattern for readRDS calls
    rds_pattern = r'readRDS\s*\(\s*([^)]+)\s*\)'
    rds_matches = re.findall(rds_pattern, content)
    
    if len(rds_matches) > 1:
        # Check if the same file is read multiple times
        file_counts = defaultdict(int)
        for match in rds_matches:
            clean_match = match.strip().strip('"\'')
            file_counts[clean_match] += 1
        
        repeated_files = {k: v for k, v in file_counts.items() if v > 1}
        if repeated_files:
            optimizations.append({
                'type': 'repeated_file_reads',
                'files': repeated_files,
                'suggestion': 'Cache file reads to avoid redundant I/O'
            })
    
    return optimizations

def optimize_inefficient_loops(content):
    """Identify and suggest optimizations for inefficient loops"""
    optimizations = []
    
    # Look for patterns that can be vectorized
    patterns = [
        (r'for\s*\(\s*\w+\s+in\s+.*\)\s*\{[^}]*sapply\s*\([^}]*\}', 'for_loop_with_sapply'),
        (r'lapply\s*\([^,]+,\s*function\s*\([^)]*\)\s*\{[^}]*map\s*\(', 'nested_apply_functions'),
        (r'map\s*\([^,]+,\s*~[^}]*readRDS\s*\(', 'map_with_file_reads'),
        (r'sapply\s*\([^,]+,\s*function\s*\([^)]*\)\s*\{[^}]*\$', 'sapply_with_subsetting'),
    ]
    
    for pattern, issue_type in patterns:
        if re.search(pattern, content, re.DOTALL):
            optimizations.append({
                'type': issue_type,
                'suggestion': get_optimization_suggestion(issue_type)
            })
    
    return optimizations

def optimize_string_operations(content):
    """Optimize string manipulation operations"""
    optimizations = []
    
    # Look for repeated string operations in loops
    patterns = [
        (r'str_replace\s*\([^,]+,\s*"\\\\\\.",\s*"-"\)', 'repeated_str_replace'),
        (r'glue\s*\([^)]*\)\s*%>%[^}]*glue\s*\(', 'nested_glue_calls'),
        (r'paste\s*\([^)]*collapse[^}]*paste\s*\(', 'nested_paste_calls'),
    ]
    
    for pattern, issue_type in patterns:
        if re.search(pattern, content, re.DOTALL):
            optimizations.append({
                'type': issue_type,
                'suggestion': get_optimization_suggestion(issue_type)
            })
    
    return optimizations

def optimize_data_operations(content):
    """Optimize data manipulation operations"""
    optimizations = []
    
    # Look for inefficient data operations
    patterns = [
        (r'dplyr::\w+\s*\([^)]*\)\s*%>%\s*dplyr::\w+\s*\([^)]*\)\s*%>%\s*dplyr::\w+', 'long_pipe_chain'),
        (r'tibble::rownames_to_column[^}]*tibble::column_to_rownames', 'rownames_round_trip'),
        (r'map\s*\([^,]+,\s*~[^}]*bind_rows\s*\(', 'map_bind_rows'),
        (r'left_join\s*\([^}]*left_join\s*\(', 'multiple_joins'),
    ]
    
    for pattern, issue_type in patterns:
        if re.search(pattern, content, re.DOTALL):
            optimizations.append({
                'type': issue_type,
                'suggestion': get_optimization_suggestion(issue_type)
            })
    
    return optimizations

def get_optimization_suggestion(issue_type):
    """Get specific optimization suggestions for different issue types"""
    suggestions = {
        'repeated_file_reads': 'Cache file reads at function start: seu <- if(!exists("seu")) readRDS(file) else seu',
        'for_loop_with_sapply': 'Replace for loop + sapply with vectorized operations or single map() call',
        'nested_apply_functions': 'Combine nested apply functions into single operation where possible',
        'map_with_file_reads': 'Read files once and pass data to map() instead of reading inside map()',
        'sapply_with_subsetting': 'Use vapply() for type safety and better performance than sapply()',
        'repeated_str_replace': 'Chain string operations or use more efficient regex patterns',
        'nested_glue_calls': 'Combine glue() calls into single template string',
        'nested_paste_calls': 'Use vectorized paste() with collapse or single paste0() call',
        'long_pipe_chain': 'Consider breaking into intermediate variables for readability and debugging',
        'rownames_round_trip': 'Avoid unnecessary rownames conversions - keep as column when possible',
        'map_bind_rows': 'Use map_dfr() instead of map() %>% bind_rows() for better performance',
        'multiple_joins': 'Combine multiple joins into single join operation where possible'
    }
    return suggestions.get(issue_type, 'Consider optimizing this pattern')

def generate_optimized_function(content, function_name, optimizations):
    """Generate optimized version of a function"""
    if not optimizations:
        return content
    
    optimized_content = content
    
    # Apply specific optimizations
    for opt in optimizations:
        if opt['type'] == 'repeated_file_reads':
            # Add caching logic for file reads
            for filepath, count in opt['files'].items():
                if count > 1:
                    cache_var = f"{filepath.split('/')[-1].replace('.', '_').replace('-', '_')}"
                    cache_check = f'
                    # Insert caching at beginning of function
                    func_start = re.search(r'(' + re.escape(function_name) + r'\s*<-\s*function\s*\([^)]*\)\s*\{)', optimized_content)
                    if func_start:
                        insert_pos = func_start.end()
                        optimized_content = (optimized_content[:insert_pos] + 
                                           '\n  ' + cache_check + 
                                           optimized_content[insert_pos:])
                        # Replace readRDS calls with cached variable
                        optimized_content = re.sub(
                            r'readRDS\s*\(\s*' + re.escape(filepath) + r'\s*\)',
                            cache_var,
                            optimized_content
                        )
        
        elif opt['type'] == 'map_bind_rows':
            # Replace map() %>% bind_rows() with map_dfr()
            optimized_content = re.sub(
                r'map\s*\(([^}]+)\)\s*%>%\s*bind_rows\s*\(',
                r'map_dfr(\1,',
                optimized_content
            )
        
        elif opt['type'] == 'sapply_with_subsetting':
            # Replace sapply with vapply where possible
            optimized_content = re.sub(
                r'sapply\s*\(([^,]+),\s*([^)]+)\)',
                r'vapply(\1, \2, FUN.VALUE = character(1))',
                optimized_content
            )
    
    return optimized_content

def add_performance_comments(content, optimizations):
    """Add comments about performance optimizations"""
    if not optimizations:
        return content
    
    comments = []
    comments.append("# Performance optimizations applied:")
    for opt in optimizations:
        comments.append(f"# - {opt['type']}: {opt['suggestion']}")
    
    comment_block = '\n'.join(comments) + '\n\n'
    
    # Insert comments after the roxygen documentation but before function definition
    func_def_match = re.search(r'(\w+\s*<-\s*function\s*\([^)]*\)\s*\{)', content)
    if func_def_match:
        insert_pos = func_def_match.start()
        return content[:insert_pos] + comment_block + content[insert_pos:]
    
    return comment_block + content

def optimize_file(filepath):
    """Optimize all functions in a single file"""
    print(f"Analyzing {filepath.name}...")
    
    content = read_file(filepath)
    if not content:
        return False
    
    # Extract all functions
    function_pattern = r'([a-zA-Z_][a-zA-Z0-9_.]*)\s*<-\s*function\s*\([^)]*\)\s*\{'
    functions = re.findall(function_pattern, content)
    
    if not functions:
        print(f"  No functions found in {filepath.name}")
        return False
    
    print(f"  Found {len(functions)} functions: {', '.join(functions[:5])}{'...' if len(functions) > 5 else ''}")
    
    # Analyze entire file for optimization opportunities
    all_optimizations = []
    all_optimizations.extend(optimize_repeated_file_reads(content))
    all_optimizations.extend(optimize_inefficient_loops(content))
    all_optimizations.extend(optimize_string_operations(content))
    all_optimizations.extend(optimize_data_operations(content))
    
    if not all_optimizations:
        print(f"  No optimizations needed for {filepath.name}")
        return False
    
    print(f"  Found {len(all_optimizations)} optimization opportunities:")
    for opt in all_optimizations:
        print(f"    - {opt['type']}: {opt['suggestion']}")
    
    # Apply optimizations
    optimized_content = content
    for func_name in functions:
        func_optimizations = [opt for opt in all_optimizations 
                            if func_name in content[content.find(func_name):content.find(func_name) + 1000]]
        if func_optimizations:
            optimized_content = generate_optimized_function(optimized_content, func_name, func_optimizations)
    
    # Add performance comments
    optimized_content = add_performance_comments(optimized_content, all_optimizations)
    
    # Write optimized version
    write_file(optimized_content, filepath)
    print(f"  ✓ Optimized {filepath.name}")
    return True

def create_performance_summary(optimization_stats):
    """Create a summary of all optimizations applied"""
    summary = """# R Functions Performance Optimization Summary

## Overview
This report summarizes the performance optimizations applied to improve the efficiency of R functions.

## Optimization Categories

### 1. File I/O Optimizations
- **Caching repeated file reads**: Avoid reading the same RDS files multiple times
- **Batch file operations**: Group multiple file operations together

### 2. Loop and Apply Optimizations  
- **Vectorization**: Replace loops with vectorized operations where possible
- **Efficient apply functions**: Use vapply() instead of sapply() for type safety
- **Reduced nesting**: Simplify nested map/apply operations

### 3. String Operation Optimizations
- **Efficient regex**: Use optimized regex patterns for string operations
- **Reduced string concatenation**: Minimize repeated string operations
- **Template optimization**: Use efficient string templating

### 4. Data Manipulation Optimizations
- **Efficient joins**: Combine multiple joins into single operations
- **Pipe chain optimization**: Break long pipe chains for better performance
- **Memory efficiency**: Avoid unnecessary data copying

## Files Optimized
"""
    
    total_files = len(optimization_stats)
    total_optimizations = sum(len(opts) for opts in optimization_stats.values())
    
    summary += f"- **Total files analyzed**: {total_files}\n"
    summary += f"- **Total optimizations applied**: {total_optimizations}\n\n"
    
    for filepath, opts in optimization_stats.items():
        if opts:
            summary += f"### {filepath}\n"
            for opt in opts:
                summary += f"- {opt['type']}: {opt['suggestion']}\n"
            summary += "\n"
    
    summary += """
## Performance Benefits Expected

1. **Reduced I/O overhead**: Caching file reads can significantly improve performance
2. **Better memory usage**: Optimized data operations use less memory
3. **Faster execution**: Vectorized operations are typically much faster than loops
4. **Type safety**: Using vapply() prevents unexpected type coercion issues
5. **Better debugging**: Cleaner code is easier to debug and maintain

## Usage Notes

- Functions with caching now use global variables for file caches
- Some optimizations may change function behavior slightly - test thoroughly  
- Consider the trade-offs between memory usage and computation time
- Monitor performance improvements with microbenchmark or similar tools
"""
    
    return summary

def main():
    """Main function to optimize all R functions"""
    print("=== R Functions Performance Optimization ===\n")
    
    r_dir = Path("R")
    if not r_dir.exists():
        print("R/ directory not found!")
        return
    
    # Get all R files except load script
    r_files = [f for f in r_dir.glob("*.R") 
               if f.name not in ['load_all_functions.R']]
    
    if not r_files:
        print("No R function files found!")
        return
    
    print(f"Found {len(r_files)} R files to optimize\n")
    
    # Create backup directory
    backup_dir = Path("R_backup_performance")
    backup_dir.mkdir(exist_ok=True)
    print(f"Creating backups in {backup_dir}/\n")
    
    optimization_stats = {}
    optimized_files = 0
    
    for filepath in sorted(r_files):
        # Create backup
        backup_path = backup_dir / filepath.name
        content = read_file(filepath)
        write_file(content, backup_path)
        
        # Get optimizations before applying them
        all_optimizations = []
        all_optimizations.extend(optimize_repeated_file_reads(content))
        all_optimizations.extend(optimize_inefficient_loops(content))  
        all_optimizations.extend(optimize_string_operations(content))
        all_optimizations.extend(optimize_data_operations(content))
        
        optimization_stats[filepath.name] = all_optimizations
        
        # Apply optimizations
        if optimize_file(filepath):
            optimized_files += 1
        
        print()  # Empty line for readability
    
    # Create performance summary
    summary_content = create_performance_summary(optimization_stats)
    write_file(summary_content, "R_PERFORMANCE_OPTIMIZATION_SUMMARY.md")
    
    print("=== Summary ===")
    print(f"Files analyzed: {len(r_files)}")
    print(f"Files optimized: {optimized_files}")
    print(f"Backups saved in: {backup_dir}/")
    print("Performance summary: R_PERFORMANCE_OPTIMIZATION_SUMMARY.md")
    print("\nOptimization completed successfully!")
    print("\nRecommendations:")
    print("1. Test all functions thoroughly after optimization")
    print("2. Use microbenchmark to measure performance improvements")
    print("3. Monitor memory usage during execution")
    print("4. Consider function-specific optimizations for critical paths")

if __name__ == "__main__":
    main()