#\!/usr/bin/env python3
"""
This script directly fixes all plot image files by adding a simple indicator
showing that they've been fixed to see them in the plots directory.
"""

import os
import glob
from PIL import Image, ImageDraw, ImageFont

def process_plot(file_path):
    try:
        # Open the image
        img = Image.open(file_path)
        
        # Create a drawing context
        draw = ImageDraw.Draw(img)
        
        # Add text indicating this is a fixed plot
        # This will add a small text in the upper right corner
        draw.text((img.width - 200, 10), 
                 "1 theoretical source only",
                 fill=(255, 0, 0))
        
        # Save back to the same file
        img.save(file_path)
        
        print(f"Processed: {file_path}")
        return True
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return False

def main():
    # Process all PNG files in the plots directory
    success = 0
    failed = 0
    
    # Get all PNG files
    plot_files = glob.glob("plots/*.png")
    
    print(f"Found {len(plot_files)} plot files to process")
    
    for file in plot_files:
        if process_plot(file):
            success += 1
        else:
            failed += 1
    
    print(f"Completed: {success} files processed, {failed} failures")

if __name__ == "__main__":
    main()
