# =============================================================================
#  dataSet.py - Consolida os labels do Roboflow (YOLO) em um unico CSV.
#
#  Cada arquivo .txt de label = um FRAME (imagem). Cada linha do arquivo = uma
#  deteccao de cacho. O frame_id identifica de qual imagem veio cada deteccao,
#  permitindo a analise espacial por frame.
#
#  Requer a pasta do export do Roboflow (YOLOv8) com as subpastas */labels/*.txt.
#  Ajuste 'base_path' para a pasta onde estao os labels.
# =============================================================================
import os
import csv

base_path = r"C:\ESD\Trabalho\Grapes.v19i.yolov8"
output_file = r"C:\ESD\Trabalho\dataset_uvas_consolidado.csv"

sample_id = 0
frame_id = 0

with open(output_file, mode="w", newline="") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(["sample_id", "frame_id", "cx", "cy", "width", "height", "area"])

    # Percorre todas as subpastas procurando arquivos de label
    for root, dirs, files in os.walk(base_path):
        for file in sorted(files):
            if file.endswith(".txt") and file not in ["classes.txt", "README.txt"]:
                txt_path = os.path.join(root, file)

                with open(txt_path, "r") as f:
                    for line in f:
                        parts = line.strip().split()
                        # YOLO: class_id cx cy w h
                        if len(parts) >= 5:
                            try:
                                cx = float(parts[1])
                                cy = float(parts[2])
                                w = float(parts[3])
                                h = float(parts[4])
                                area = w * h
                                writer.writerow([sample_id, frame_id, cx, cy, w, h, area])
                                sample_id += 1
                            except ValueError:
                                pass

                frame_id += 1  # cada arquivo de label e um frame

print(f"Sucesso! {sample_id} deteccoes em {frame_id} frames -> {output_file}")
