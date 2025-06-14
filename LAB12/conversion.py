from moviepy.editor import VideoFileClip
from pathlib import Path

# Folder z plikami mp4
input_folder = Path("plots")

# Pobierz wszystkie pliki .mp4 w folderze
mp4_files = list(input_folder.glob("*.mp4"))

print(f"Znaleziono {len(mp4_files)} plików .mp4")

for mp4_file in mp4_files:
    gif_path = mp4_file.with_suffix(".gif")
    print(f"Konwertuję: {mp4_file.name} → {gif_path.name}")

    # Wczytaj i konwertuj do GIF
    clip = VideoFileClip(str(mp4_file))
    clip.write_gif(str(gif_path), fps=10)
    print(f"Zapisano: {gif_path}\n")
