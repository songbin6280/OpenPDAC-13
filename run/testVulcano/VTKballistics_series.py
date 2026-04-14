import os
import pandas as pd
import json 

def write_vtk(points, diameters, velocities, output_file):
    n_points = len(points)
    with open(output_file, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Ballistics data\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write(f"POINTS {n_points} float\n")
        for p in points:
            f.write(f"{p[0]} {p[1]} {p[2]}\n")

        f.write(f"\nPOINT_DATA {n_points}\n")
        f.write("SCALARS diameter float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for d in diameters:
            f.write(f"{d}\n")

        f.write("VECTORS velocity float\n")
        for v in velocities:
            f.write(f"{v[0]} {v[1]} {v[2]}\n")


def process_simulation(sim_path):
    """
    Processa una simulazione: legge i file CSV in postProcessing/cloudInfo1/<time>
    e genera i corrispondenti VTK direttamente in postProcessing/cloudInfo1
    """
    cloudinfo_path = os.path.join(sim_path, "postProcessing", "cloudInfo1")
    print(f"🔎 Cerco in: {cloudinfo_path}")
    if not os.path.exists(cloudinfo_path):
        print(f"⚠️ Nessuna cartella cloudInfo1 in {sim_path}")
        return

    # Lista per raccogliere i dati per il file .series
    series_files = []

    # Ordina le cartelle numericamente
    folders = sorted(os.listdir(cloudinfo_path), key=lambda x: float(x) if x.replace('.','',1).isdigit() else 0)
    for folder in folders:
        folder_path = os.path.join(cloudinfo_path, folder)
        if not os.path.isdir(folder_path):
            continue

        csv_path = os.path.join(folder_path, "output.csv")
        if not os.path.exists(csv_path):
            continue

        df = pd.read_csv(csv_path)

        required_cols = {"x", "y", "z", "d", "Ux", "Uy", "Uz"}
        if not required_cols.issubset(df.columns):
            print(f"⚠️ Mancano colonne in {csv_path}")
            continue

        points = df[["x", "y", "z"]].values
        diameters = df["d"].values
        velocities = df[["Ux", "Uy", "Uz"]].values

        # Nome del file VTK
        vtk_filename = f"particles_{folder}.vtk"
        vtk_path = os.path.join(cloudinfo_path, vtk_filename)
        write_vtk(points, diameters, velocities, vtk_path)
        print(f"✅ Scritto: {vtk_path}")

        # Aggiunta alla lista per il file .series
        try:
            time_val = float(folder)
            series_files.append({
                "name": vtk_filename,
                "time": time_val
                })
        except ValueError:
            continue
    # Scrittura del file particles.vtk.series
    if series_files:
        series_path = os.path.join(cloudinfo_path, "particles.vtk.series")
        series_data = {
                "file-series-version": "1.0",
                "files": series_files
                }
        with open(series_path, "w") as f:
            json.dump(series_data, f, indent=4)
        print(f"\n🚀 File serie creato: {series_path}")
        print("Ora in ParaView puoi aprire direttamente 'particles.vtk.series' per vedere l'animazione.")



if __name__ == "__main__":
    base_path = os.getcwd() #usa la cartella corrente
    process_simulation(base_path)



