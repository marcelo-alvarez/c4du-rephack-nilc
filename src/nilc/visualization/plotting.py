"""Plotting utilities for NILC maps."""
import numpy as np
import os
from PIL import Image, ImageDraw, ImageFont
from pixell import enmap, enplot
from nilc.io.map_io import load_act_map, get_stokes_component


def plot_benchmark_map(benchmark_path, output_path=None, downsample_factor=20):
    """
    Visualize a benchmark CMB temperature map with proper labels and formatting.
    
    Parameters
    ----------
    benchmark_path : str
        Path to the benchmark FITS file
    output_path : str, optional
        Output path for the plot. If None, saves to data/output/benchmark_visualization.png
    downsample_factor : int, default=20
        Factor by which to downsample the map for visualization
        
    Returns
    -------
    str
        Path to the saved plot file
    """
    print(f"Loading benchmark map: {benchmark_path}")
    map_data = load_act_map(benchmark_path)
    temp_map = get_stokes_component(map_data, 'I')
    
    print(f"Map shape: {temp_map.shape}")
    print(f"Map WCS: {temp_map.wcs}")
    print(f"Statistics: min={np.nanmin(temp_map):.3e}, max={np.nanmax(temp_map):.3e}, "
          f"mean={np.nanmean(temp_map):.3e}, std={np.nanstd(temp_map):.3e}")
    
    # Downsample for visualization
    temp_ds = temp_map[::downsample_factor, ::downsample_factor]
    # Preserve WCS information by creating new enmap
    temp_ds_enmap = enmap.ndmap(temp_ds, temp_map.wcs)
    print(f"Downsampled shape: {temp_ds_enmap.shape}")
    
    # Set output path
    if output_path is None:
        output_dir = "data/output"
        os.makedirs(output_dir, exist_ok=True)
        output_prefix = os.path.join(output_dir, "benchmark_visualization")
    else:
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        output_prefix = os.path.splitext(output_path)[0]
    
    # Use pixell enplot to create visualization with proper CAR projection
    # Enable ticks and labels (enplot will use WCS to generate coordinate labels)
    plot_obj = enplot.plot(temp_ds_enmap, colorbar=True, downgrade=2, ticks=2)
    
    # Set white background and add labels by modifying the PIL image
    for plot in plot_obj:
        if hasattr(plot, 'img') and plot.img:
            # Convert RGBA to RGB with white background
            img = plot.img
            if img.mode == 'RGBA':
                # Create white background
                background = Image.new('RGB', img.size, (255, 255, 255))
                # Composite the image onto white background
                background.paste(img, mask=img.split()[3])  # Use alpha channel as mask
                plot.img = background
            elif img.mode in ['RGB', 'L']:
                # Already RGB or grayscale, just ensure it's RGB
                if img.mode == 'L':
                    plot.img = img.convert('RGB')
            
            # Add padding around the image for labels
            padding_top = 60
            padding_bottom = 50
            padding_left = 80
            padding_right = 80
            
            img_width, img_height = plot.img.size
            
            # Estimate plot width vs colorbar width (enplot typically puts colorbar on right)
            # We'll assume the colorbar is roughly 10% of the width, so plot is ~90%
            # But to be safe, we'll add extra space after the main plot area
            # For now, let's add the padding and then place labels correctly
            plot_area_width = int(img_width * 0.85)  # Estimate plot takes ~85% of width
            colorbar_start = plot_area_width
            
            new_width = img_width + padding_left + padding_right
            new_height = img_height + padding_top + padding_bottom
            
            # Create new image with white background and padding
            padded_img = Image.new('RGB', (new_width, new_height), (255, 255, 255))
            # Paste original image in the center
            padded_img.paste(plot.img, (padding_left, padding_top))
            plot.img = padded_img
            
            # Now add labels using PIL ImageDraw
            draw = ImageDraw.Draw(plot.img)
            
            # Try to use a default font, fallback to basic if not available
            try:
                font_medium = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans.ttf", 14)
            except:
                try:
                    font_medium = ImageFont.load_default()
                except:
                    font_medium = None
            
            # X-axis label (bottom center, directly under plot, not under colorbar) - Right Ascension
            x_label = "Right Ascension (degrees)"
            # Center under the plot area (not the entire image including colorbar)
            # Plot center is at: padding_left + plot_area_width/2
            plot_center_x = padding_left + plot_area_width // 2
            if font_medium:
                bbox = draw.textbbox((0, 0), x_label, font=font_medium)
                text_width = bbox[2] - bbox[0]
                draw.text((plot_center_x - text_width // 2, padding_top + img_height + 10), x_label, 
                         fill=(0, 0, 0), font=font_medium)
            else:
                draw.text((plot_center_x - 100, padding_top + img_height + 10), x_label, fill=(0, 0, 0))
            
            # Y-axis label (left side, rotated 90 degrees)
            y_label = "Declination (degrees)"
            if font_medium:
                # Create a temporary image for rotated text
                temp_img = Image.new('RGBA', (200, 200), (255, 255, 255, 0))
                temp_draw = ImageDraw.Draw(temp_img)
                bbox = temp_draw.textbbox((0, 0), y_label, font=font_medium)
                text_width = bbox[2] - bbox[0]
                text_height = bbox[3] - bbox[1]
                temp_draw.text((0, 0), y_label, fill=(0, 0, 0, 255), font=font_medium)
                # Rotate 90 degrees counterclockwise
                rotated_text = temp_img.rotate(90, expand=True, fillcolor=(255, 255, 255, 0))
                # Paste rotated text on left side
                rot_width, rot_height = rotated_text.size
                plot.img.paste(rotated_text, (padding_left // 2 - rot_width // 2, 
                                             padding_top + img_height // 2 - rot_height // 2), 
                              rotated_text)
            else:
                draw.text((10, padding_top + img_height // 2 - 50), y_label, fill=(0, 0, 0))
            
            # Colorbar label (right side, rotated 90 degrees) - Temperature with units
            cbar_label = "Temperature (μK)"
            if font_medium:
                # Create a temporary image for rotated text with enough space
                temp_img = Image.new('RGBA', (300, 300), (255, 255, 255, 0))
                temp_draw = ImageDraw.Draw(temp_img)
                bbox = temp_draw.textbbox((0, 0), cbar_label, font=font_medium)
                text_width = bbox[2] - bbox[0]
                text_height = bbox[3] - bbox[1]
                # Center text in the temporary image
                temp_draw.text((150 - text_width // 2, 150 - text_height // 2), cbar_label, 
                              fill=(0, 0, 0, 255), font=font_medium)
                # Rotate 90 degrees counterclockwise
                rotated_text = temp_img.rotate(90, expand=True, fillcolor=(255, 255, 255, 0))
                # Paste rotated text on right side near colorbar (after the plot area)
                rot_width, rot_height = rotated_text.size
                # Position it to the right of the plot, near where colorbar would be
                cbar_label_x = padding_left + plot_area_width + 30  # Space after plot, near colorbar
                plot.img.paste(rotated_text, (cbar_label_x - rot_width // 2, 
                                             padding_top + img_height // 2 - rot_height // 2), 
                              rotated_text)
            else:
                cbar_label_x = padding_left + plot_area_width + 30
                draw.text((cbar_label_x, padding_top + img_height // 2 - 50), cbar_label, fill=(0, 0, 0))
    
    enplot.write(output_prefix, plot_obj)
    # enplot adds .png extension, so the final file is output_prefix + ".png"
    final_output_path = output_prefix + ".png"
    print(f"Saved plot to: {final_output_path}")
    
    return final_output_path


if __name__ == "__main__":
    # Benchmark path from validation instructions
    benchmark_path = "/secret/path/to/target/output/target_T.fits"
    plot_benchmark_map(benchmark_path)

