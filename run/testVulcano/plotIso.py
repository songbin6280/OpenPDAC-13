#
# =============================================================================
#
#          PyVista Video Generator for OpenFOAM Simulations
#
# Description:
#   This script automates the creation of high-quality MP4 videos from
#   OpenFOAM simulation outputs. It visualizes the time evolution of
#   isosurfaces and Lagrangian particles over a static 3D terrain.
#
# Features:
#   - Renders a static terrain from a .vtk file.
#   - Applies a texture, stretched to fit the surface perfectly, from a .jpg file.
#   - Renders two evolving isosurfaces from .vtk file sequences.
#   - Renders evolving Lagrangian particles (ballistics) from .csv files,
#     with particle size represented by a colormap.
#   - Uses multiprocessing to significantly speed up frame generation.
#   - Optional interactive mode to set up the base camera view and texture.
#   - Command-line arguments to control framerate, parallel processes,
#     particle size scaling, color mapping, and colorbar visibility.
#   - Automatically assembles rendered frames into MP4 videos using FFMPEG.
#
# =============================================================================
#

import pyvista as pv
import glob
import os
import re
import time
import numpy as np
import argparse
import shutil
import subprocess
from multiprocessing import Pool, cpu_count
import pandas as pd

# ----- DEFAULT SETTINGS -----

# Default frame rate for the output videos (frames per second).
DEFAULT_FRAMERATE = 10
# Default resolution for the output videos (width, height) in pixels.
DEFAULT_VIDEO_RESOLUTION = (1920, 1080)
# Default number of parallel processes to use for rendering.
# Uses all CPU cores minus one to leave resources for the system.
DEFAULT_NUM_PROCESSES = max(1, cpu_count() - 1)
# Background color for the rendering window.
BACKGROUND_COLOR = 'white'
# Default colormap for particle sizing. 'viridis' is perceptually uniform.
DEFAULT_PARTICLE_CMAP = 'viridis'

# ----- ARGUMENT PARSER -----
# Sets up the command-line interface to configure the script's behavior.
parser = argparse.ArgumentParser(
    description="""
Generate videos from OpenFOAM surface and Lagrangian particle outputs.

This tool searches for simulation results in the 'postProcessing' directory
and creates several video sequences from different camera angles, including
static isometric views, a top-down view, and a 360-degree orbit animation.
""",
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(
    "--interactive",
    action="store_true",
    help=("Enable interactive view to rotate texture, save camera "
          "(for isometric base), and then proceed."),
)
parser.add_argument(
    "--np", "--num_processes",
    type=int,
    default=None,
    help=(f"Number of parallel processes for frame generation "
          f"(default: {DEFAULT_NUM_PROCESSES} or auto-detected)."),
)
parser.add_argument(
    "--fr", "--framerate",
    type=int,
    default=None,
    help=f"Framerate for the output videos (default: {DEFAULT_FRAMERATE}).")

parser.add_argument(
    "--ps", "--particle_scale",
    type=float,
    default=1.0,
    help=("Scaling factor for ballistic particle diameters "
          "(e.g., 10 to make them 10x larger). Default: 1.0")
)
parser.add_argument(
    "--cmap",
    type=str,
    default=DEFAULT_PARTICLE_CMAP,
    help=(f"Colormap for particle diameters (e.g., viridis, jet, coolwarm). "
          f"Default: {DEFAULT_PARTICLE_CMAP}")
)
parser.add_argument(
    "--ncb", "--no_colorbar",
    action="store_true",
    dest='no_colorbar', # Explicitly define the destination variable
    help="If set, disables the particle size colorbar in the videos."
)
args = parser.parse_args()


# ----- SETTINGS FROM ARGS OR DEFAULTS -----

# Framerate for the final videos.
FRAMERATE = args.fr if args.fr is not None else DEFAULT_FRAMERATE
# Number of parallel processes for rendering. Capped by the number of CPU cores.
NUM_PROCESSES = args.np if args.np is not None else DEFAULT_NUM_PROCESSES
NUM_PROCESSES = max(1, min(NUM_PROCESSES, cpu_count()))
# Resolution of the final videos.
VIDEO_RESOLUTION = DEFAULT_VIDEO_RESOLUTION
# Multiplicative factor for particle diameters.
PARTICLE_SCALE_FACTOR = args.ps
# Colormap for particle diameters.
PARTICLE_CMAP = args.cmap
# Flag for colorbar visibility.
SHOW_COLORBAR = not args.no_colorbar

# ----- GLOBAL CAMERA AND CONTROL VARIABLES -----

# Flag to signal continuation after the interactive session.
continue_animation_flag = False
# Stores the camera position set by the user in interactive mode.
user_defined_isometric_base_view = None
# Stores the calculated top-down camera position.
calculated_top_down_view = None

# ----- HELPER FUNCTIONS -----

def read_foam_dict_value(filepath, keyword):
    """
    Reads a numerical value from a simplified OpenFOAM dictionary file.
    """
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        content = re.sub(r'//.*', '', content)
        content = re.sub(r'/\*.*?\*/', '', content, flags=re.DOTALL)
        match = re.search(rf'{keyword}\s+([\d.eE+-]+);', content)
        if match:
            return float(match.group(1))
        return None
    except Exception:
        return None

def key_press_callback_interactive(plotter_instance):
    """
    Callback function executed when a key is pressed in the interactive plotter.
    """
    global continue_animation_flag, user_defined_isometric_base_view
    user_defined_isometric_base_view = plotter_instance.camera_position
    print(f"User-defined isometric base view saved: {user_defined_isometric_base_view}")
    continue_animation_flag = True
    plotter_instance.exit()

def rotate_view_around_focal_point(camera_view, angle_deg):
    """
    Calculates a new camera position by rotating around the focal point.
    """
    pos, focal, view_up = camera_view
    angle_rad = np.deg2rad(angle_deg)
    vec_focal_to_pos = np.array(pos) - np.array(focal)
    R_z = np.array([[np.cos(angle_rad), -np.sin(angle_rad), 0],
                    [np.sin(angle_rad), np.cos(angle_rad), 0],
                    [0, 0, 1]])
    vec_focal_to_pos_rotated = R_z @ vec_focal_to_pos
    new_pos = np.array(focal) + vec_focal_to_pos_rotated
    new_view_up = R_z @ np.array(view_up)
    return (new_pos.tolist(), focal, new_view_up.tolist())

def normalize_texture_coordinates(mesh):
    """
    Normalizes the texture coordinates of a mesh to the exact [0, 1] range.
    """
    if "Texture Coordinates" not in mesh.point_data:
        return mesh
    t_coords = mesh.point_data['Texture Coordinates']
    u_min, u_max = t_coords[:, 0].min(), t_coords[:, 0].max()
    if (u_max - u_min) != 0:
        t_coords[:, 0] = (t_coords[:, 0] - u_min) / (u_max - u_min)
    v_min, v_max = t_coords[:, 1].min(), t_coords[:, 1].max()
    if (v_max - v_min) != 0:
        t_coords[:, 1] = (t_coords[:, 1] - v_min) / (v_max - v_min)
    mesh.point_data['Texture Coordinates'] = t_coords
    print("Texture coordinates have been normalized to prevent repeating.")
    return mesh

# ----- PART 1: TERRAIN/TEXTURE LOADING & CAMERA SETUP -----

TERRAIN_VTK_PATH = "postProcessing/terrain/0/terrain.vtk"
TEXTURE_IMAGE_PATH = "constant/DEM/texture.jpg"
print("Loading terrain mesh...")
if not os.path.exists(TERRAIN_VTK_PATH):
    raise FileNotFoundError(f"FATAL: Terrain file not found: {TERRAIN_VTK_PATH}")
terrain_mesh = pv.read(TERRAIN_VTK_PATH)
if not terrain_mesh or terrain_mesh.n_points == 0:
    raise ValueError(f"FATAL: Terrain mesh from '{TERRAIN_VTK_PATH}' is empty or invalid.")
print(f"Terrain mesh loaded. Bounds: {terrain_mesh.bounds}, Center: {terrain_mesh.center}")

terrain_mesh.compute_normals(auto_orient_normals=True, consistent_normals=True, inplace=True)
terrain_mesh.texture_map_to_plane(inplace=True)
terrain_mesh = normalize_texture_coordinates(terrain_mesh)

texture_exists_globally = os.path.exists(TEXTURE_IMAGE_PATH)
terrain_color_fallback = "lightgray"
if texture_exists_globally:
    print(f"Texture file '{TEXTURE_IMAGE_PATH}' found.")
else:
    print(f"Warning: Texture file '{TEXTURE_IMAGE_PATH}' not found. Using solid color.")

texture_rotation_angle_deg = [0]
texture_flip_x = [False]
texture_flip_y = [False]

def apply_texture_transformations_to_coords(mesh_with_uvs):
    """
    Applies interactive rotation and flipping to the texture coordinates.
    """
    if "Texture Coordinates" not in mesh_with_uvs.point_data:
        mesh_with_uvs.texture_map_to_plane(inplace=True)
        mesh_with_uvs = normalize_texture_coordinates(mesh_with_uvs)
    t_coords = mesh_with_uvs.point_data['Texture Coordinates'].copy()
    t_coords[:, 0] -= 0.5; t_coords[:, 1] -= 0.5
    angle_rad = np.deg2rad(texture_rotation_angle_deg[0])
    u, v = t_coords[:, 0].copy(), t_coords[:, 1].copy()
    t_coords[:, 0] = np.cos(angle_rad) * u - np.sin(angle_rad) * v
    t_coords[:, 1] = np.sin(angle_rad) * u + np.cos(angle_rad) * v
    t_coords[:, 0] += 0.5; t_coords[:, 1] += 0.5
    if texture_flip_x[0]: t_coords[:, 0] = 1.0 - t_coords[:, 0]
    if texture_flip_y[0]: t_coords[:, 1] = 1.0 - t_coords[:, 1]
    mesh_with_uvs.point_data['Texture Coordinates'] = t_coords
    return mesh_with_uvs

if texture_exists_globally:
    terrain_mesh = apply_texture_transformations_to_coords(terrain_mesh)

def get_specific_view(mesh, view_type='iso'):
    """
    Calculates a standard camera position for a given mesh.
    """
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(mesh.copy())
    if view_type == 'iso':
        plotter.camera_position = 'iso'; plotter.reset_camera(render=False)
    elif view_type == 'xy':
        plotter.view_xy(negative=False)
    cam_pos = plotter.camera_position
    plotter.close(); del plotter
    return cam_pos

if args.interactive:
    print("\nStarting interactive setup for ISOMETRIC BASE VIEW...")
    if texture_exists_globally:
        interactive_texture_object = pv.read_texture(TEXTURE_IMAGE_PATH)
    plotter_interactive = pv.Plotter(window_size=VIDEO_RESOLUTION, off_screen=False)
    plotter_interactive.background_color = BACKGROUND_COLOR
    if texture_exists_globally and interactive_texture_object:
        actor_terrain_interactive = plotter_interactive.add_mesh(terrain_mesh, texture=interactive_texture_object, smooth_shading=True)
    else:
        actor_terrain_interactive = plotter_interactive.add_mesh(terrain_mesh, color=terrain_color_fallback, smooth_shading=True)
    plotter_interactive.add_axes(zlabel='Z (Up)', xlabel='X', ylabel='Y')
    initial_iso_cam = get_specific_view(terrain_mesh, 'iso')
    plotter_interactive.camera_position = initial_iso_cam
    def cb_rotate_texture_interactive():
        texture_rotation_angle_deg[0] = (texture_rotation_angle_deg[0] + 90) % 360
        global terrain_mesh
        terrain_mesh_copy = terrain_mesh.copy(); terrain_mesh_copy.texture_map_to_plane(inplace=True)
        terrain_mesh_copy = normalize_texture_coordinates(terrain_mesh_copy)
        terrain_mesh = apply_texture_transformations_to_coords(terrain_mesh_copy)
        plotter_interactive.update_scalars(terrain_mesh.point_data['Texture Coordinates'], mesh=actor_terrain_interactive, render=True)
    def cb_flip_texture_x_interactive():
        texture_flip_x[0] = not texture_flip_x[0]
        global terrain_mesh
        terrain_mesh_copy = terrain_mesh.copy(); terrain_mesh_copy.texture_map_to_plane(inplace=True)
        terrain_mesh_copy = normalize_texture_coordinates(terrain_mesh_copy)
        terrain_mesh = apply_texture_transformations_to_coords(terrain_mesh_copy)
        plotter_interactive.update_scalars(terrain_mesh.point_data['Texture Coordinates'], mesh=actor_terrain_interactive, render=True)
    def cb_flip_texture_y_interactive():
        texture_flip_y[0] = not texture_flip_y[0]
        global terrain_mesh
        terrain_mesh_copy = terrain_mesh.copy(); terrain_mesh_copy.texture_map_to_plane(inplace=True)
        terrain_mesh_copy = normalize_texture_coordinates(terrain_mesh_copy)
        terrain_mesh = apply_texture_transformations_to_coords(terrain_mesh_copy)
        plotter_interactive.update_scalars(terrain_mesh.point_data['Texture Coordinates'], mesh=actor_terrain_interactive, render=True)
    plotter_interactive.add_key_event("r", cb_rotate_texture_interactive)
    plotter_interactive.add_key_event("x", cb_flip_texture_x_interactive)
    plotter_interactive.add_key_event("y", cb_flip_texture_y_interactive)
    plotter_interactive.add_key_event("Return", lambda: key_press_callback_interactive(plotter_interactive))
    print("  Interactive Window: Adjust ISOMETRIC base view. R:Rotate Txt, X/Y:Flip Txt, Enter:Save Cam, Q:Quit.")
    plotter_interactive.show(title="Interactive Isometric Base View Setup")
    if not continue_animation_flag:
        user_defined_isometric_base_view = get_specific_view(terrain_mesh, 'iso')
    if plotter_interactive.renderer:
        plotter_interactive.close()
    del plotter_interactive
    if interactive_texture_object:
        del interactive_texture_object
else:
    print("\nNon-interactive mode. Using default calculated isometric base view.")
    user_defined_isometric_base_view = get_specific_view(terrain_mesh, 'iso')

print("Calculating top-down view...")
calculated_top_down_view = get_specific_view(terrain_mesh, 'xy')
print(f"--- User-defined/Default ISOMETRIC BASE view for animations: {user_defined_isometric_base_view} ---")
print(f"--- Calculated TOP-DOWN view for animations: {calculated_top_down_view} ---")

# ----- PART 2: PREPARING ANIMATION DATA -----

print("\nLocating VTK and CSV files for animation...")
ISO_SURFACES_DIR_PATTERN = "postProcessing/isoSurfaces/*"
PARTICLES_DIR_PATTERN = "postProcessing/cloudInfo1/*"
ISO4_VTK_NAME = "iso4.vtk"; ISO6_VTK_NAME = "iso6.vtk"; PARTICLE_CSV_NAME = "output.csv"

def get_timestep_from_dir_path(dir_path):
    try: return float(os.path.basename(dir_path))
    except (ValueError, TypeError): return -1

timestep_dirs = sorted(glob.glob(ISO_SURFACES_DIR_PATTERN), key=get_timestep_from_dir_path)
valid_timestep_dirs = [d for d in timestep_dirs if get_timestep_from_dir_path(d) >= 0]
iso4_vtk_files, iso6_vtk_files, simulation_timesteps = [], [], []
for ts_dir in valid_timestep_dirs:
    iso4_path, iso6_path = os.path.join(ts_dir, ISO4_VTK_NAME), os.path.join(ts_dir, ISO6_VTK_NAME)
    if os.path.exists(iso4_path) and os.path.exists(iso6_path):
        iso4_vtk_files.append(iso4_path); iso6_vtk_files.append(iso6_path)
        simulation_timesteps.append(get_timestep_from_dir_path(ts_dir))
if not iso4_vtk_files: raise FileNotFoundError(f"No complete VTK sets found in {ISO_SURFACES_DIR_PATTERN}.")
print(f"Found {len(iso4_vtk_files)} timesteps with both iso4 and iso6 VTK files.")

print("Locating particle CSV files...")
particle_dirs = glob.glob(PARTICLES_DIR_PATTERN)
particle_data_map = {}
for p_dir in particle_dirs:
    timestep = get_timestep_from_dir_path(p_dir)
    if timestep >= 0:
        csv_path = os.path.join(p_dir, PARTICLE_CSV_NAME)
        if os.path.exists(csv_path): particle_data_map[timestep] = csv_path

if not particle_data_map:
    print("Warning: No particle CSV files found.")
    global_particle_diameter_range = None
else:
    print(f"Found particle data for {len(particle_data_map)} timesteps.")
    print("Calculating global particle diameter range for consistent colormap...")
    global_min_d, global_max_d = float('inf'), float('-inf')
    for csv_file in particle_data_map.values():
        try:
            df = pd.read_csv(csv_file)
            if not df.empty and 'd' in df.columns:
                global_min_d = min(global_min_d, df['d'].min()); global_max_d = max(global_max_d, df['d'].max())
        except Exception as e: print(f"  - Warning: Could not process {csv_file}: {e}")
    if global_min_d > global_max_d:
        global_particle_diameter_range = None
        print("  - Could not determine a valid diameter range.")
    else:
        global_particle_diameter_range = (global_min_d, global_max_d)
        print(f"  - Global diameter range: {global_particle_diameter_range[0]:.4f} to {global_particle_diameter_range[1]:.4f}")

OUTPUT_FRAMES_ROOT_DIR = "postProcessing/frames_png"; OUTPUT_VIDEOS_ROOT_DIR = "postProcessing/videos_mp4"
os.makedirs(OUTPUT_FRAMES_ROOT_DIR, exist_ok=True); os.makedirs(OUTPUT_VIDEOS_ROOT_DIR, exist_ok=True)

# ----- FUNCTION FOR RENDERING A SINGLE FRAME (WORKER) -----

def render_single_frame_worker(frame_index, num_total_frames, iso4_filepath,
                               iso6_filepath, current_sim_time,
                               camera_view_for_frame, output_png_filepath,
                               passed_terrain_mesh, passed_texture_exists_flag,
                               passed_texture_image_path,
                               passed_terrain_color,
                               particle_csv_filepath,
                               particle_scale_factor,
                               particle_cmap,
                               particle_clim,
                               show_colorbar_flag):
    """
    Renders a single frame of the animation.
    """
    frame_plotter = pv.Plotter(off_screen=True, window_size=VIDEO_RESOLUTION)
    frame_plotter.background_color = BACKGROUND_COLOR
    if passed_texture_exists_flag:
        try:
            worker_texture_object = pv.read_texture(passed_texture_image_path)
            frame_plotter.add_mesh(passed_terrain_mesh, texture=worker_texture_object, smooth_shading=True, name="terrain")
        except Exception as e:
            if frame_index == 0: print(f"[WORKER {os.getpid()}] Error loading texture: {e}")
            frame_plotter.add_mesh(passed_terrain_mesh, color=passed_terrain_color, smooth_shading=True, name="terrain")
    elif passed_terrain_mesh:
        frame_plotter.add_mesh(passed_terrain_mesh, color=passed_terrain_color, smooth_shading=True, name="terrain")
    try:
        mesh4, mesh6 = pv.read(iso4_filepath), pv.read(iso6_filepath)
        if mesh4.n_points > 0: frame_plotter.add_mesh(mesh4, color="tan", opacity=0.7, smooth_shading=True, name="iso4")
        if mesh6.n_points > 0: frame_plotter.add_mesh(mesh6, color="lightblue", opacity=0.7, smooth_shading=True, name="iso6")
    except Exception as e:
        if frame_index == 0: print(f"[WORKER {os.getpid()}] Error loading VTK: {e}")

    if particle_csv_filepath and os.path.exists(particle_csv_filepath):
        try:
            df = pd.read_csv(particle_csv_filepath)
            if not df.empty and all(col in df.columns for col in ['x', 'y', 'z', 'd']):
                points, diameters = df[['x', 'y', 'z']].values, df['d'].values
                particle_cloud = pv.PolyData(points)
                particle_cloud['diameter_color'] = diameters
                particle_cloud['diameter_glyph'] = diameters * particle_scale_factor
                sphere_glyph = pv.Sphere(theta_resolution=8, phi_resolution=8)
                glyphs = particle_cloud.glyph(scale='diameter_glyph', geom=sphere_glyph, orient=False)
                
                # --- CORRECTED SECTION ---
                
                # 1. Add the mesh and get the single ACTOR object in return
                actor = frame_plotter.add_mesh(
                    glyphs,
                    scalars='diameter_color',
                    cmap=particle_cmap,
                    clim=particle_clim,
                    show_scalar_bar=False,
                    name="particles"
                )

                # 2. If requested, add the scalar bar manually, getting the mapper FROM THE ACTOR
                if show_colorbar_flag and particle_clim is not None:
                    title = "Particle Diameter (m)"
                    if particle_scale_factor != 1.0:
                        title += f" [Visual Size x{particle_scale_factor:.0f}]"
                                            
                    sbar_width = 0.5
                    frame_plotter.add_scalar_bar(
                        title=title,
                        mapper=actor.mapper,  # <-- The crucial fix is here
                        n_labels=3,
                        italic=False,
                        bold=False,
                        title_font_size=20,
                        label_font_size=18,
                        color='black',
                        width=sbar_width,
                        position_x=0.5 - (sbar_width / 2),
                        position_y=0.05
                    )
                # --- END OF CORRECTED SECTION ---
                
        except Exception as e:
            if frame_index == 0: print(f"[WORKER {os.getpid()}] Error processing particles: {e}")

    frame_plotter.add_text(f"Time: {current_sim_time:.2f}s", font_size=15, position="upper_right", color="black", shadow=True)
    frame_plotter.add_axes(zlabel='Z', xlabel='X', ylabel='Y')
    frame_plotter.camera_position = camera_view_for_frame
    try:
        frame_plotter.screenshot(output_png_filepath, transparent_background=False)
    except Exception as e:
        print(f"[WORKER {os.getpid()}] EXCEPTION during screenshot: {e}")
    frame_plotter.close(); del frame_plotter
    return output_png_filepath

# ----- FUNCTION TO ASSEMBLE VIDEO -----

def assemble_video_from_frames(frames_input_dir, video_output_path, framerate_val, delete_frames_after=True):
    """
    Assembles a sequence of PNG frames into an MP4 video using FFMPEG.
    """
    print(f"Assembling video: '{video_output_path}' from '{frames_input_dir}' at {framerate_val} FPS.")
    png_filename_pattern = os.path.join(frames_input_dir, "frame_%04d.png")
    ffmpeg_command = ['ffmpeg', '-r', str(framerate_val), '-i', png_filename_pattern, '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-vf', "pad=ceil(iw/2)*2:ceil(ih/2)*2", '-an', '-y', video_output_path]
    try:
        subprocess.run(ffmpeg_command, check=True, capture_output=True, text=True)
        print(f"Video successfully saved: {video_output_path}")
        if delete_frames_after:
            shutil.rmtree(frames_input_dir)
            print(f"Directory '{frames_input_dir}' removed.")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: FFMPEG failed.\n  Command: {' '.join(e.cmd)}\n  Stderr:\n{e.stderr}")
        print(f"PNG frames preserved in '{frames_input_dir}'.")
    except FileNotFoundError:
        print("ERROR: FFMPEG not found. Ensure it's installed and in your PATH.")
        print(f"PNG frames preserved in '{frames_input_dir}'.")

# ----- MAIN FUNCTION TO GENERATE A VIDEO SEQUENCE -----

def create_animation_sequence(base_camera_view, rotation_angle_degrees=None, use_circular_orbit=False, sequence_base_name="video_seq"):
    """
    Manages the entire process of rendering and assembling a video sequence.
    """
    current_sequence_frames_dir = os.path.join(OUTPUT_FRAMES_ROOT_DIR, sequence_base_name)
    os.makedirs(current_sequence_frames_dir, exist_ok=True)
    print(f"\nPreparing frames for sequence: '{sequence_base_name}' in '{current_sequence_frames_dir}'")
    num_animation_frames = len(iso4_vtk_files)
    if num_animation_frames == 0:
        print("No frames to render. Skipping."); return
    worker_tasks_list = []
    for i in range(num_animation_frames):
        current_sim_time = simulation_timesteps[i]
        current_camera_for_frame = base_camera_view
        if use_circular_orbit:
            orbit_angle_deg = 360.0 * i / num_animation_frames
            current_camera_for_frame = rotate_view_around_focal_point(base_camera_view, orbit_angle_deg)
        elif rotation_angle_degrees is not None:
            current_camera_for_frame = rotate_view_around_focal_point(base_camera_view, rotation_angle_degrees)
        output_png_filepath = os.path.join(current_sequence_frames_dir, f"frame_{i:04d}.png")
        particle_file_for_frame = particle_data_map.get(current_sim_time, None)
        task_arg_tuple = (
            i, num_animation_frames, iso4_vtk_files[i], iso6_vtk_files[i],
            current_sim_time, current_camera_for_frame, output_png_filepath,
            terrain_mesh, texture_exists_globally, TEXTURE_IMAGE_PATH,
            terrain_color_fallback if not texture_exists_globally else None,
            particle_file_for_frame, PARTICLE_SCALE_FACTOR, PARTICLE_CMAP,
            global_particle_diameter_range, SHOW_COLORBAR
        )
        worker_tasks_list.append(task_arg_tuple)
    print(f"Rendering {num_animation_frames} frames for '{sequence_base_name}' using {NUM_PROCESSES} processes...")
    start_time = time.time()
    with Pool(processes=NUM_PROCESSES) as frame_pool:
        results = frame_pool.starmap(render_single_frame_worker, worker_tasks_list)
    end_time = time.time()
    print(f"Frames for '{sequence_base_name}' rendered in {end_time - start_time:.2f}s.")
    rendered_files_count = sum(1 for r in results if r and os.path.exists(r))
    if rendered_files_count != num_animation_frames:
        print(f"Warning: Expected {num_animation_frames} frames, got {rendered_files_count}.")
    final_video_path = os.path.join(OUTPUT_VIDEOS_ROOT_DIR, f"{sequence_base_name}.mp4")
    assemble_video_from_frames(current_sequence_frames_dir, final_video_path, FRAMERATE)

# ----- MAIN EXECUTION SCRIPT -----

if __name__ == "__main__":
    print(f"Script started. Using {NUM_PROCESSES} processes for rendering at {FRAMERATE} FPS.")
    print(f"Particle diameter scale factor: {PARTICLE_SCALE_FACTOR}x")
    print(f"Particle colormap: '{PARTICLE_CMAP}'")
    print(f"Show particle colorbar: {SHOW_COLORBAR}")
    print(f"Video resolution: {VIDEO_RESOLUTION[0]}x{VIDEO_RESOLUTION[1]}, Background: {BACKGROUND_COLOR}")
    if not user_defined_isometric_base_view or not calculated_top_down_view:
        raise RuntimeError("FATAL: Base camera views were not properly initialized.")

    create_animation_sequence(base_camera_view=calculated_top_down_view, sequence_base_name="video_view_top_down")
    static_isometric_angles = [0, 90, 180, 270]
    for angle in static_isometric_angles:
        create_animation_sequence(base_camera_view=user_defined_isometric_base_view, rotation_angle_degrees=angle, sequence_base_name=f"video_view_iso_{angle}deg")
    create_animation_sequence(base_camera_view=user_defined_isometric_base_view, use_circular_orbit=True, sequence_base_name="video_orbit_iso_360deg")

    print("\nAll requested video sequences processed.")
    print(f"Final videos should be in: '{OUTPUT_VIDEOS_ROOT_DIR}'")
