# Proof of concept while developing ideas for this pipeline
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from typing import List
from shiny import App, ui, render, reactive, req

# ---------- Data Normalization ----------
def read_and_normalize_files(file_paths: List[str], label_prefix: str) -> pd.DataFrame:
    merged_df = None

    for path in file_paths:
        ext = os.path.splitext(path)[1].lower()
        delimiter = "\t" if ext == ".tsv" else ","

        try:
            df = pd.read_csv(path, header=None, delimiter=delimiter, usecols=[0, 1])
        except Exception as e:
            print(f"⚠️ Could not read {path}: {e}")
            continue

        df = df.dropna()
        try:
            df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')
            df = df.dropna()
        except:
            continue

        if df.shape[1] < 2:
            continue

        norm_col = (df.iloc[:, 1] - df.iloc[:, 1].min()) / (df.iloc[:, 1].max() - df.iloc[:, 1].min())
        label = f"{label_prefix}_{os.path.splitext(os.path.basename(path))[0]}"

        temp_df = pd.DataFrame({
            "repeat length": df.iloc[:, 0],
            label: norm_col
        })

        if merged_df is None:
            merged_df = temp_df
        else:
            merged_df = pd.merge(merged_df, temp_df, on="repeat length", how="inner")

    if merged_df is None:
        raise ValueError(f"No valid files found in {label_prefix}")

    return merged_df

# ---------- File Resolver ----------
def resolve_uploaded_paths(fileinfo_list):
    if not fileinfo_list:
        return []

    paths = []
    for item in fileinfo_list:
        path = item["datapath"]
        if os.path.isdir(path):
            files = glob.glob(os.path.join(path, "*.histogram")) + \
                    glob.glob(os.path.join(path, "*.csv")) + \
                    glob.glob(os.path.join(path, "*.tsv")) + \
                    glob.glob(os.path.join(path, "*.txt"))
            paths.extend(files)
        else:
            paths.append(path)
    return paths

# ---------- Plot Generator ----------
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from textwrap import wrap

def generate_delta_plot(dataframe, output_path="delta_plot.pdf", bin_size=2, title="Delta Plot", color="red"):
    df = dataframe.copy()
    df = df.sort_values("repeat length")
    df_filtered = df.iloc[::bin_size, :]

    # Set A4 portrait dimensions in inches (8.27 × 11.69)
    fig = plt.figure(figsize=(8.27, 11.69))
    gs = GridSpec(2, 1, height_ratios=[1, 1], figure=fig)

    # --- Plot Area ---
    ax1 = fig.add_subplot(gs[0])
    ax1.axhline(0, color='gray', linestyle='--', linewidth=1)
    ax1.plot(df_filtered["repeat length"], df_filtered["delta"], linestyle='-', marker='o',
             markersize=4, color=color)
    ax1.set_title(title, fontsize=14)
    ax1.set_xlabel("Repeat Length", fontsize=12)
    ax1.set_ylabel("Delta (Day42_avg - Day0_avg)", fontsize=12)
    ax1.set_xlim(0, 200)

    # --- Description Area (bottom of page, not an Axes) ---
    description_title = "Description of Plot:"
    description_body = (
        "This plot shows the average change in CAG repeat size between Day 0 and Day 42. "
        "The X-axis represents the CAG repeat length bins, and the Y-axis represents the delta "
        "(Day 42 average minus Day 0 average). "
        "An increased delta indicates a greater decrease in CAG repeat size over time, suggesting contraction of the repeat region. "
        "Each point is based on normalized input histograms."
    )

    # Wrap the body text to fit the page
    wrapped_body = "\n".join(wrap(description_body, width=100))

    # Add title and body as fig-level text
    fig.text(0.05, 0.42, description_title, fontsize=12, weight="bold", ha="left")
    fig.text(0.05, 0.39, wrapped_body, fontsize=10, ha="left", va="top")

    # Save and close
    plt.tight_layout(rect=[0, 0.45, 1, 1])  # Leave room at bottom for description
    plt.savefig(output_path)
    plt.close()


# ---------- UI ----------
app_ui = ui.page_fluid(
    ui.h2("Normalize and Compare Histograms"),
    ui.input_file("day0_input", "Upload Day 0 Files or Folder", multiple=True, accept=[".csv", ".tsv", ".txt", ".histogram"]),
    ui.input_file("day42_input", "Upload Day 42 Files or Folder", multiple=True, accept=[".csv", ".tsv", ".txt", ".histogram"]),
    ui.input_slider("bin_size", "Bin Size (Point Interval)", min=1, max=20, value=2),
    ui.input_text("plot_title", "Plot Title", placeholder="Enter plot title...", value="Delta Plot"),
    ui.input_select("plot_color", "Plot Color", choices=["red", "blue", "green", "black", "orange", "pink", "purple", "grey"], selected="red"),
    ui.input_action_button("normalise_button", "Normalise", class_="btn-success", style="margin-top: 20px;"),
    ui.input_action_button("plot_button", "Generate Delta Plot", class_="btn-warning", style="margin-top: 20px; margin-left: 10px;"),
    ui.output_ui("status_text")
)



# ---------- Server ----------
def server(input, output, session):
    merged_data = reactive.Value(None)

    @reactive.effect
    @reactive.event(input.normalise_button)
    def _():
        day0_files = resolve_uploaded_paths(input.day0_input())
        day42_files = resolve_uploaded_paths(input.day42_input())

        if not day0_files or not day42_files:
            print("⚠️ One or both file groups missing.")
            return

        df_day0 = read_and_normalize_files(day0_files, "Day0")
        df_day42 = read_and_normalize_files(day42_files, "Day42")

        combined = pd.merge(df_day0, df_day42, on="repeat length", how="inner")

        day0_cols = [col for col in combined.columns if col.startswith("Day0_")]
        day42_cols = [col for col in combined.columns if col.startswith("Day42_")]

        combined["Day0_avg"] = combined[day0_cols].mean(axis=1)
        combined["Day42_avg"] = combined[day42_cols].mean(axis=1)
        combined["delta"] = combined["Day42_avg"] - combined["Day0_avg"]

        cols = [col for col in combined.columns if col != "delta"]
        combined = combined[cols + ["delta"]]

        combined.to_csv("normalized_with_delta.csv", index=False)
        merged_data.set(combined)
        print("✅ CSV saved to: normalized_with_delta.csv")

    @reactive.effect
    @reactive.event(input.plot_button)
    def _():
        df = merged_data.get()
        if df is not None:
            bin_size = input.bin_size()
            title = input.plot_title()
            color = input.plot_color()
            generate_delta_plot(df, output_path="delta_plot.pdf", bin_size=bin_size, title=title, color=color)
            print(f"✅ Delta plot saved to: delta_plot.pdf")

    @output
    @render.ui
    def status_text():
        if merged_data.get() is not None:
            return ui.p("✅ Normalized data and delta plot saved to normalised_with_deslta.csv")
        return ui.p("📂 Upload files and click 'Normalise'.")
    

# ---------- Run App ----------
app = App(app_ui, server)
