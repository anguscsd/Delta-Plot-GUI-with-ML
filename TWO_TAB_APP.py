import wx
import os
import glob
import platform
import subprocess
import threading
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from textwrap import wrap
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages
from pubsub import pub

# ---------Create "created by" banner ----------------
def created_by_banner(parent, text="Created by Ruban Rex and Angus Dixon"):
    banner_panel = wx.Panel(parent) 
    banner_panel.SetBackgroundColour('#D3D3D3')  # banner background colour
    sizer = wx.BoxSizer(wx.HORIZONTAL)

    credit_text = wx.StaticText(banner_panel, label=text)
    credit_text.SetForegroundColour('#333333')  # text colour
    credit_text.SetFont(wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL))

    sizer.AddStretchSpacer()
    sizer.Add(credit_text, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
    sizer.AddStretchSpacer()

    banner_panel.SetSizer(sizer)
    return banner_panel

# Shared helpers
# ------------------------
def resolve_uploaded_paths(paths):
    resolved = []
    for path in paths:
        if os.path.isdir(path):
            files = glob.glob(os.path.join(path, "*.histogram")) + \
                    glob.glob(os.path.join(path, "*.csv")) + \
                    glob.glob(os.path.join(path, "*.tsv")) + \
                    glob.glob(os.path.join(path, "*.txt"))
            resolved.extend(files)
        else:
            resolved.append(path)
    return resolved

# ---------- Analyzer: read & normalize ----------
def read_and_normalize_files(file_paths, label_prefix):
    merged_df = None

    for path in file_paths:
        ext = os.path.splitext(path)[1].lower()
        is_hist = (ext == ".histogram")

        # Choose a default delimiter
        # .histogram and .tsv are tab-separated; others default to comma
        delimiter = "\t" if ext in (".tsv", ".histogram") else ","

        # .histogram files have a 6-line header
        skip = 6 if is_hist else 0

        # First attempt: read with our expected delimiter
        try:
            df = pd.read_csv(
                path,
                header=None,
                usecols=[0, 1],
                sep=delimiter,
                skiprows=skip,
                engine="python",
            )
        except Exception:
            # Fallback: let pandas sniff the delimiter
            try:
                df = pd.read_csv(
                    path,
                    header=None,
                    usecols=[0, 1],
                    sep=None,           # automatic delimiter detection
                    skiprows=skip,
                    engine="python",
                )
            except Exception as e2:
                print(f"⚠️ Could not read {path}: {e2}")
                continue

        # Clean & coerce
        df = df.dropna()
        df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors="coerce")
        df = df.dropna()

        if df.shape[1] < 2 or df.empty:
            continue

        # Normalise with divide-by-zero guard
        y = df.iloc[:, 1]
        denom = y.max() - y.min()
        if not np.isfinite(denom) or denom == 0:
            continue

        norm_col = (y - y.min()) / denom
        label = f"{label_prefix}_{os.path.splitext(os.path.basename(path))[0]}"

        temp_df = pd.DataFrame({
            "repeat length": df.iloc[:, 0],
            label: norm_col
        })

        if merged_df is None:
            merged_df = temp_df
        else:
            merged_df = pd.merge(merged_df, temp_df, on="repeat length", how="inner")

    if merged_df is None or merged_df.empty:
        raise ValueError(f"No valid files found in {label_prefix} (check delimiter/format).")

    return merged_df


# ---------- Analyzer: PDF plot (supports error bars) ----------
def generate_delta_plot(dataframe, output_path="delta_plot.pdf", bin_size=2, title="Delta Plot",
                        color="red", error_mode="None", dataframe2=None, color2="black"):
    df = dataframe.copy().sort_values("repeat length")
    df_f = df.iloc[::bin_size, :]

    with PdfPages(output_path) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        gs = GridSpec(2, 1, height_ratios=[1, 1], figure=fig)
        ax1 = fig.add_subplot(gs[0])

        ax1.axhline(0, color='gray', linestyle='--', linewidth=1)

        # --- Set 1 ---
        x1 = df_f["repeat length"]
        y1 = df_f["delta"]
        yerr1 = None
        if error_mode != "None":
            d0c = [c for c in dataframe.columns if c.startswith("Day0_")]
            d42c = [c for c in dataframe.columns if c.startswith("Day42_")]
            if d0c and d42c:
                d0v = dataframe[d0c].iloc[::bin_size, :]
                d42v = dataframe[d42c].iloc[::bin_size, :]
                if error_mode == "Standard Deviation (SD)":
                    yerr1 = np.sqrt(d0v.std(axis=1)**2 + d42v.std(axis=1)**2)
                else:
                    yerr1 = np.sqrt(d0v.sem(axis=1)**2 + d42v.sem(axis=1)**2)

        if yerr1 is not None:
            ax1.errorbar(x1, y1, yerr=yerr1, fmt='o-', color=color, capsize=3, markersize=3, label="Set 1")
        else:
            ax1.plot(x1, y1, linestyle='-', marker='o', markersize=4, color=color, label="Set 1")

        # --- Optional Set 2 ---
        if dataframe2 is not None:
            df2 = dataframe2.copy().sort_values("repeat length")
            df2_f = df2.iloc[::bin_size, :]
            x2 = df2_f["repeat length"]
            y2 = df2_f["delta"]

            yerr2 = None
            if error_mode != "None":
                d0c2 = [c for c in dataframe2.columns if c.startswith("Day0_")]
                d42c2 = [c for c in dataframe2.columns if c.startswith("Day42_")]
                if d0c2 and d42c2:
                    d0v2 = dataframe2[d0c2].iloc[::bin_size, :]
                    d42v2 = dataframe2[d42c2].iloc[::bin_size, :]
                    if error_mode == "Standard Deviation (SD)":
                        yerr2 = np.sqrt(d0v2.std(axis=1)**2 + d42v2.std(axis=1)**2)
                    else:
                        yerr2 = np.sqrt(d0v2.sem(axis=1)**2 + d42v2.sem(axis=1)**2)

            if yerr2 is not None:
                ax1.errorbar(x2, y2, yerr=yerr2, fmt='s--', color=color2, capsize=3, markersize=3, label="Set 2")
            else:
                ax1.plot(x2, y2, linestyle='--', marker='s', markersize=4, color=color2, label="Set 2")

        ax1.set_title(title, fontsize=14)
        ax1.set_xlabel("Repeat Length", fontsize=12)
        ax1.set_ylabel("Delta (Day42_avg - Day0_avg)", fontsize=12)
        ax1.set_xlim(30, 200)
        ax1.legend()

        # Description / footer
        description_title = "Description of Plot:"
        description_body = (
            "The delta plot shows the average change in repeat size between Day 0 and Day 42. "
            "If two datasets are provided, both are overlaid to facilitate comparison."
        )
        wrapped_body = "\n".join(wrap(description_body, width=100))
        fig.text(0.05, 0.42, description_title, fontsize=12, weight="bold", ha="left")
        fig.text(0.05, 0.39, wrapped_body, fontsize=10, ha="left", va="top")

        plt.tight_layout(rect=[0, 0.45, 1, 1])
        pdf.savefig(fig)
        plt.close(fig)

        # Optional page comparing Day0/Day42 averages (overlay both sets)
        fig2, ax2 = plt.subplots(figsize=(8.27, 4))
        x = df["repeat length"]
        if "Day0_avg" in df and "Day42_avg" in df:
            ax2.plot(x, df["Day0_avg"], label="Day 0 (Set 1)")
            ax2.plot(x, df["Day42_avg"], label="Day 42 (Set 1)")
        if dataframe2 is not None and "Day0_avg" in dataframe2 and "Day42_avg" in dataframe2:
            x2 = dataframe2["repeat length"]
            ax2.plot(x2, dataframe2["Day0_avg"], linestyle="--", label="Day 0 (Set 2)")
            ax2.plot(x2, dataframe2["Day42_avg"], linestyle="--", label="Day 42 (Set 2)")
        ax2.set_xlabel("Repeat Length")
        ax2.set_ylabel("Normalized Frequency")
        ax2.set_title("Histogram Comparison: Day 0 vs Day 42")
        ax2.set_xlim(0, 200)
        ax2.legend()
        plt.tight_layout()
        pdf.savefig(fig2)
        plt.close(fig2)

        # Quick Day0/Day42 average overlay
        hist_fig = plot_histogram_comparison(dataframe)
        pdf.savefig(hist_fig)
        plt.close(hist_fig)

def plot_histogram_comparison(df):
    fig, ax = plt.subplots(figsize=(8.27, 4))
    x = df["repeat length"]
    if "Day0_avg" in df and "Day42_avg" in df:
        y_day0 = df["Day0_avg"]
        y_day42 = df["Day42_avg"]
        ax.plot(x, y_day0, label="Day 0", linewidth=2)
        ax.plot(x, y_day42, label="Day 42", linewidth=2)
        ax.legend()
    ax.set_xlabel("Repeat Length")
    ax.set_ylabel("Normalized Frequency")
    ax.set_title("Histogram Comparison: Day 0 vs Day 42")
    ax.set_xlim(0, 200)
    plt.tight_layout()
    return fig

# -------------------------
# TAB 1: Generate via Docker
# -------------------------
class GenerateHistogramsPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # Add credit banner at top
        vbox.Add(created_by_banner(self), flag=wx.EXPAND | wx.TOP | wx.BOTTOM, border=5)

        # Docker .tar file dir (as in your original)
        row_tar = wx.BoxSizer(wx.HORIZONTAL)
        row_tar.Add(wx.StaticText(self, label='Docker .tar File Directory:'), 0, wx.RIGHT|wx.ALIGN_CENTER_VERTICAL, 8)
        self.tarDirCtrl = wx.TextCtrl(self)
        row_tar.Add(self.tarDirCtrl, 1, wx.EXPAND)
        tarBrowseBtn = wx.Button(self, label='Browse', size=(80, 30))
        tarBrowseBtn.Bind(wx.EVT_BUTTON, self.on_browse_tar)
        row_tar.Add(tarBrowseBtn, 0, wx.LEFT, 8)
        vbox.Add(row_tar, 0, wx.EXPAND | wx.ALL, 8)

        # FASTA files directory
        row_fa = wx.BoxSizer(wx.HORIZONTAL)
        row_fa.Add(wx.StaticText(self, label='FASTA Files Directory:'), 0, wx.RIGHT|wx.ALIGN_CENTER_VERTICAL, 8)
        self.fastaDirCtrl = wx.TextCtrl(self)
        row_fa.Add(self.fastaDirCtrl, 1, wx.EXPAND)
        faBrowse = wx.Button(self, label='Browse', size=(80, 30))
        faBrowse.Bind(wx.EVT_BUTTON, self.on_browse_fasta)
        row_fa.Add(faBrowse, 0, wx.LEFT, 8)
        vbox.Add(row_fa, 0, wx.EXPAND | wx.ALL, 8)

        # Output base directory (NEW)
        row_out = wx.BoxSizer(wx.HORIZONTAL)
        row_out.Add(wx.StaticText(self, label='Histogram Output Base Folder:'), 0, wx.RIGHT|wx.ALIGN_CENTER_VERTICAL, 8)
        self.outDirCtrl = wx.TextCtrl(self)
        row_out.Add(self.outDirCtrl, 1, wx.EXPAND)
        outBrowse = wx.Button(self, label='Browse', size=(80, 30))
        outBrowse.Bind(wx.EVT_BUTTON, self.on_browse_out)
        row_out.Add(outBrowse, 0, wx.LEFT, 8)
        vbox.Add(row_out, 0, wx.EXPAND | wx.ALL, 8)

        # Profile and optional instability index
        row_prof = wx.BoxSizer(wx.HORIZONTAL)
        row_prof.Add(wx.StaticText(self, label='Profile Mode:'), 0, wx.RIGHT|wx.ALIGN_CENTER_VERTICAL, 8)
        self.profile_choice = wx.Choice(self, choices=["restrictive", "permissive"])
        self.profile_choice.SetSelection(0)
        row_prof.Add(self.profile_choice, 0)
        self.instability_checkbox = wx.CheckBox(self, label='Run Instability Index')
        row_prof.Add(self.instability_checkbox, 0, wx.LEFT, 12)
        vbox.Add(row_prof, 0, wx.ALL, 8)

        # Run
        self.run_btn = wx.Button(self, label='Run RD Program', size=(160, 38))
        self.run_btn.Bind(wx.EVT_BUTTON, self.on_run)
        vbox.Add(self.run_btn, 0, wx.ALL, 8)

        self.gauge = wx.Gauge(self, range=100)
        vbox.Add(self.gauge, 0, wx.EXPAND | wx.ALL, 8)

        vbox.Add(wx.StaticText(self, label="Run log:"), 0, wx.LEFT|wx.TOP, 8)
        self.outputCtrl = wx.TextCtrl(self, style=wx.TE_MULTILINE | wx.TE_READONLY)
        vbox.Add(self.outputCtrl, 1, wx.EXPAND | wx.ALL, 8)

        self.SetSizer(vbox)

    def on_browse_tar(self, _):
        with wx.DirDialog(self, "Choose the Docker .tar file directory:") as d:
            if d.ShowModal() == wx.ID_OK:
                self.tarDirCtrl.SetValue(d.GetPath())

    def on_browse_fasta(self, _):
        with wx.DirDialog(self, "Choose the FASTA files directory:") as d:
            if d.ShowModal() == wx.ID_OK:
                self.fastaDirCtrl.SetValue(d.GetPath())

    def on_browse_out(self, _):
        with wx.DirDialog(self, "Choose the output base directory:") as d:
            if d.ShowModal() == wx.ID_OK:
                self.outDirCtrl.SetValue(d.GetPath())

    def log(self, msg):
        wx.CallAfter(self.outputCtrl.AppendText, msg + "\n")

    def _classify_day_mainthread(self, filename):
        """Run on the main thread: show a dialog to map file → Day0/Day42 if ambiguous."""
        name = os.path.basename(filename).lower()
        if "day42" in name or "day_42" in name or "d42" in name:
            return "Day42"
        if "day0" in name or "day_0" in name or "d0" in name:
            return "Day0"
        dlg = wx.SingleChoiceDialog(self, f"Assign sample '{os.path.basename(filename)}' to:",
                                    "Choose Group", ["Day0", "Day42"])
        try:
            if dlg.ShowModal() == wx.ID_OK:
                return dlg.GetStringSelection()
        finally:
            dlg.Destroy()
        return "Day0"


    def on_run(self, _):
        tar_dir = self.tarDirCtrl.GetValue()
        fasta_dir = self.fastaDirCtrl.GetValue()
        out_base = self.outDirCtrl.GetValue()
        profile = self.profile_choice.GetStringSelection()
        run_instability = self.instability_checkbox.GetValue()

        if not tar_dir or not os.path.isdir(tar_dir):
            wx.MessageBox("Invalid Docker .tar file directory.", "Error", wx.OK|wx.ICON_ERROR); return
        if not fasta_dir or not os.path.isdir(fasta_dir):
            wx.MessageBox("Invalid FASTA files directory.", "Error", wx.OK|wx.ICON_ERROR); return
        if not out_base:
            wx.MessageBox("Please choose an output base directory.", "Error", wx.OK|wx.ICON_ERROR); return

        # Create Day0/Day42 subfolders
        day0_dir = os.path.join(out_base, "Day0_histograms")
        day42_dir = os.path.join(out_base, "Day42_histograms")
        os.makedirs(day0_dir, exist_ok=True)
        os.makedirs(day42_dir, exist_ok=True)

        # PRE-CLASSIFY FASTAs ON THE MAIN THREAD
        fastas = [f for f in os.listdir(fasta_dir) if f.lower().endswith(('.fasta', '.fa'))]
        mapping = {}
        for fasta_file in fastas:
            mapping[fasta_file] = self._classify_day_mainthread(fasta_file)  # <-- new method below

        self.run_btn.Disable()
        self.gauge.Pulse()
        self.outputCtrl.Clear()
        self.log("Starting Docker runs...")

        # Pass mapping to worker
        threading.Thread(
            target=self._worker,
            args=(fasta_dir, day0_dir, day42_dir, profile, run_instability, mapping),
            daemon=True
        ).start()


    def _worker(self, fasta_dir, day0_dir, day42_dir, profile, run_instability, mapping):
        try:
            prf_file = ("/app/RepeatDetector/Profiles/CAG/Annex10_cag_correctedFreq_notlog_AND_Complete.prf"
                        if profile == "restrictive"
                        else "/app/RepeatDetector/Profiles/CAG/Annex2_cag.prf")

            fastas = [f for f in os.listdir(fasta_dir) if f.lower().endswith(('.fasta', '.fa'))]
            total = len(fastas)
            produced_day0, produced_day42 = [], []

            for i, fasta_file in enumerate(fastas, 1):
                group = mapping.get(fasta_file, "Day0")  # Use pre-chosen group
                out_dir = day0_dir if group == "Day0" else day42_dir

                file_name = os.path.splitext(fasta_file)[0]
                output_basename = f"{file_name}.rest" if profile == "restrictive" else f"{file_name}.per"

                docker_cmd = f"""/usr/local/bin/docker run --rm \
-v "{fasta_dir}:/mnt/fasta" \
-v "{out_dir}:/mnt/out" \
-e LD_LIBRARY_PATH=/app/RepeatDetector/build/external/htslib/src/htslib \
repeat-detector \
/usr/local/bin/RepeatDetecter \
--prf {prf_file} /mnt/fasta/{fasta_file} \
--output-name "/mnt/out/{output_basename}" -o histogram --with-revcomp --cycle-range [0:500] --verbose"""

                self.log(f"$ {docker_cmd}")
                p = subprocess.Popen(docker_cmd, shell=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True)
                for line in p.stdout:
                    self.log(line.rstrip())
                rc = p.wait()
                if rc != 0:
                    self.log(f"❌ Error: Docker returned {rc} for {fasta_file}")
                    continue

                hist_path = os.path.join(out_dir, f"{output_basename}.histogram")
                if os.path.isfile(hist_path):
                    (produced_day0 if group == "Day0" else produced_day42).append(hist_path)
                    if run_instability:
                        try:
                            self._run_instability_index(hist_path)
                        except Exception as e:
                            self.log(f"Instability index error: {e}")

                wx.CallAfter(self.gauge.SetValue, int(100 * i / max(1, total)))

            self.log(f"✅ Finished. Day0 files: {len(produced_day0)}, Day42 files: {len(produced_day42)}")
            # Publish to analysis tab
            from pubsub import pub  # if you’ve switched to pypubsub
            wx.CallAfter(pub.sendMessage, "histograms_ready",
                        day0=produced_day0, day42=produced_day42,
                        day0_dir=day0_dir, day42_dir=day42_dir)

        except Exception as e:
            self.log(f"❌ Fatal error: {e}")
        finally:
            wx.CallAfter(self.gauge.SetValue, 0)
            wx.CallAfter(self.run_btn.Enable)

    def _run_instability_index(self, histogram_path):
        try:
            warnings.simplefilter(action='ignore', category=FutureWarning)

            if not os.path.isfile(histogram_path):
                self.log(f"Histogram file not found: {histogram_path}")
                return

            # Load histogram data
            df = pd.read_table(histogram_path, skiprows=range(0, 6), names=["repeat_size", "reads"])
            sample_name = os.path.basename(histogram_path).replace('.rest.histogram', '').replace('.per.histogram', '')

            mode = df.loc[df['reads'].idxmax(), 'repeat_size']
            highest_peak = df['reads'].max()
            threshold = 0.05 * highest_peak
            sub_df = df[df['reads'] > threshold].copy()

            sub_df['normalized_peak_height'] = sub_df['reads'] / sub_df['reads'].sum()
            max_index = sub_df['reads'].idxmax()
            sub_df['change_from_main_allele'] = sub_df.index - max_index
            sub_df['normalized_peak'] = sub_df['normalized_peak_height'] * sub_df['change_from_main_allele']
            instability_index = sub_df['normalized_peak'].sum()

            contraction_index = sub_df[sub_df['change_from_main_allele'] < 0]['normalized_peak'].sum()
            expansion_index = sub_df[sub_df['change_from_main_allele'] > 0]['normalized_peak'].sum()

            # Thread-safe, non-GUI plotting using Agg (no pyplot)
            fig = Figure(figsize=(10, 6))
            ax = fig.add_subplot(111)
            ax.bar(df['repeat_size'], df['reads'], edgecolor='black', label="Reads")
            ax.axhline(y=threshold, linestyle='--', label="Threshold")
            ax.set_title(f"{sample_name}: Instability Index")
            ax.set_xlabel("Repeat Size")
            ax.set_ylabel("Reads")
            ax.legend()

            plot_path = f"{os.path.splitext(histogram_path)[0]}_plot.png"
            canvas = FigureCanvasAgg(fig)
            fig.savefig(plot_path)
            # (no plt.close needed)

            # Write results
            output_dir = os.path.dirname(histogram_path)
            output_txt_path = os.path.join(output_dir, "results.txt")
            with open(output_txt_path, mode='a') as txtfile:
                txtfile.write(f"Sample: {sample_name}\n")
                txtfile.write(f"Mode: {mode}\n")
                txtfile.write(f"Instability Index: {round(instability_index, 2)}\n")
                txtfile.write(f"Contraction Index: {round(contraction_index, 2)}\n")
                txtfile.write(f"Expansion Index: {round(expansion_index, 2)}\n")
                txtfile.write(f"__________\n")

            # Log (marshal to UI thread already via self.log)
            self.log(f"Sample: {sample_name}")
            self.log(f"Mode: {mode}")
            self.log(f"Instability Index: {round(instability_index, 2)}")
            self.log(f"Contraction Index: {round(contraction_index, 2)}")
            self.log(f"Expansion Index: {round(expansion_index, 2)}")
            self.log(f"Plot saved as {plot_path}")
            self.log(f"Results saved to {output_txt_path}")

        except Exception as e:
            self.log(f"Error running instability index: {e}")


# -------------------------
# TAB 2: Analyze (your existing GUI refactored into a Panel)
# -------------------------
class AnalyzeHistogramsPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)

        vbox = wx.BoxSizer(wx.VERTICAL)

        # Add credit banner at top
        vbox.Add(created_by_banner(self), flag=wx.EXPAND | wx.TOP | wx.BOTTOM, border=5)

        self.day0_btn = wx.Button(self, label="Select Day 0 Files/Folders")
        self.day42_btn = wx.Button(self, label="Select Day 42 Files/Folders")
        self.day0_paths = []
        self.day42_paths = []
        
        # --- (Set 2, optional) ---
        self.day0_btn2 = wx.Button(self, label="Select Day 0 Files/Folders (Set 2)")
        self.day42_btn2 = wx.Button(self, label="Select Day 42 Files/Folders (Set 2)")

        self.day0_paths2 = []
        self.day42_paths2 = []

        self.day0_status2 = wx.StaticText(self, label="📁 No Day 0 files selected (Set 2).")
        self.day42_status2 = wx.StaticText(self, label="📁 No Day 42 files selected (Set 2).")

        self.title_txt = wx.TextCtrl(self, value="Delta Plot", size=(400, -1))
        self.color_choice = wx.Choice(self, choices=["red", "blue", "green", "black", "orange", "pink", "purple", "grey"])
        self.color_choice.SetSelection(0)

        self.error_bar_choice = wx.Choice(self, choices=["None", "Standard Deviation (SD)", "Standard Error (SEM)"])
        self.error_bar_choice.SetSelection(0)

        self.bin_spin = wx.SpinCtrl(self, min=1, max=20, initial=2)

        self.norm_btn = wx.Button(self, label="Normalise")
        self.plot_btn = wx.Button(self, label="Generate Delta Plot")

        self.status = wx.StaticText(self, label="📂 Upload files and click 'Normalise'.")
        self.day0_status = wx.StaticText(self, label="📁 No Day 0 files selected.")
        self.day42_status = wx.StaticText(self, label="📁 No Day 42 files selected.")

        self.csv_path_ctrl = wx.TextCtrl(self, value="normalized_with_delta.csv", size=(400, -1))
        self.csv_browse_btn = wx.Button(self, label="Browse CSV Output")

        self.pdf_path_ctrl = wx.TextCtrl(self, value="delta_plot.pdf", size=(400, -1))
        self.pdf_browse_btn = wx.Button(self, label="Browse PDF Output")

        self.csv_browse_btn.Bind(wx.EVT_BUTTON, self.on_browse_csv)
        self.pdf_browse_btn.Bind(wx.EVT_BUTTON, self.on_browse_pdf)

        self.view_plot_btn = wx.Button(self, label="View Delta Plot")
        self.view_plot_btn.Bind(wx.EVT_BUTTON, self.on_view_plot)

        self.plot_panel = wx.Panel(self, size=(550, 300))
        self.plot_canvas = None

        # Layout (same as your working frame)
        vbox.Add(wx.StaticText(self, label="Plot Title:"), 0, wx.LEFT | wx.TOP, 10)
        vbox.Add(self.title_txt, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 10)

        vbox.Add(wx.StaticText(self, label="Plot Color:"), 0, wx.LEFT | wx.TOP, 10)
        vbox.Add(self.color_choice, 0, wx.LEFT, 10)

        vbox.Add(wx.StaticText(self, label="Error Bars (Optional):"), 0, wx.LEFT | wx.TOP, 10)
        vbox.Add(self.error_bar_choice, 0, wx.LEFT, 10)

        vbox.Add(wx.StaticText(self, label="Bin Size (Point Interval):"), 0, wx.LEFT | wx.TOP, 10)
        vbox.Add(self.bin_spin, 0, wx.LEFT, 10)

        vbox.Add(wx.StaticText(self, label="CSV Output File:"), 0, wx.LEFT | wx.TOP, 10)
        h_csv = wx.BoxSizer(wx.HORIZONTAL)
        h_csv.Add(self.csv_path_ctrl, 1, wx.LEFT | wx.RIGHT, 5)
        h_csv.Add(self.csv_browse_btn, 0, wx.RIGHT, 10)
        vbox.Add(h_csv, 0, wx.EXPAND)

        vbox.Add(wx.StaticText(self, label="PDF Plot Output File:"), 0, wx.LEFT | wx.TOP, 10)
        h_pdf = wx.BoxSizer(wx.HORIZONTAL)
        h_pdf.Add(self.pdf_path_ctrl, 1, wx.LEFT | wx.RIGHT, 5)
        h_pdf.Add(self.pdf_browse_btn, 0, wx.RIGHT, 10)
        vbox.Add(h_pdf, 0, wx.EXPAND)

        vbox.Add(self.day0_btn, 0, wx.EXPAND | wx.ALL, 10)
        vbox.Add(self.day0_status, 0, wx.LEFT | wx.RIGHT, 10)

        vbox.Add(self.day42_btn, 0, wx.EXPAND | wx.ALL, 10)
        vbox.Add(self.day42_status, 0, wx.LEFT | wx.RIGHT, 10)

        vbox.Add(wx.StaticLine(self), 0, wx.EXPAND | wx.ALL, 5)
        vbox.Add(wx.StaticText(self, label="Optional: Second Dataset"), 0, wx.LEFT | wx.TOP, 10)

        vbox.Add(self.day0_btn2, 0, wx.EXPAND | wx.ALL, 10)
        vbox.Add(self.day0_status2, 0, wx.LEFT | wx.RIGHT, 10)

        vbox.Add(self.day42_btn2, 0, wx.EXPAND | wx.ALL, 10)
        vbox.Add(self.day42_status2, 0, wx.LEFT | wx.RIGHT, 10)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.norm_btn, 0, wx.RIGHT, 5)
        hbox.Add(self.plot_btn, 0)
        vbox.Add(hbox, 0, wx.LEFT | wx.BOTTOM, 10)

        vbox.Add(self.status, 0, wx.ALL, 10)

        vbox.Add(self.view_plot_btn, 0, wx.LEFT | wx.BOTTOM, 10)
        vbox.Add(wx.StaticLine(self), 0, wx.EXPAND | wx.ALL, 5)
        vbox.Add(wx.StaticText(self, label="Embedded Plot Viewer:"), 0, wx.LEFT, 10)
        vbox.Add(self.plot_panel, 1, wx.EXPAND | wx.ALL, 10)

        self.SetSizer(vbox)

        # Bindings
        self.day0_btn.Bind(wx.EVT_BUTTON, self.on_day0_select)
        self.day42_btn.Bind(wx.EVT_BUTTON, self.on_day42_select)
        self.norm_btn.Bind(wx.EVT_BUTTON, self.on_normalize)
        self.plot_btn.Bind(wx.EVT_BUTTON, self.on_plot)
        self.day0_btn2.Bind(wx.EVT_BUTTON, self.on_day0_select2)
        self.day42_btn2.Bind(wx.EVT_BUTTON, self.on_day42_select2)

        self.merged_df = None
        self.merged_df2 = None

        # Subscribe to outputs from Tab 1
        pub.subscribe(self._on_histos_ready, "histograms_ready")

    # ----- Tab-to-tab handoff -----
    def _on_histos_ready(self, day0, day42, day0_dir, day42_dir):
        # Accept lists of histogram file paths from Tab 1
        self.day0_paths = day0[:] if day0 else [day0_dir]
        self.day42_paths = day42[:] if day42 else [day42_dir]
        self.day0_status.SetLabel(f"📥 From Docker: {len(day0)} Day 0 files in {os.path.basename(day0_dir)}")
        self.day42_status.SetLabel(f"📥 From Docker: {len(day42)} Day 42 files in {os.path.basename(day42_dir)}")
        self.status.SetLabel("Ready to normalise. Click 'Normalise'.")

        # Auto-switch to this tab
        nb = self.GetParent()
        if isinstance(nb, wx.Notebook):
            nb.SetSelection(1)

    # ----- File pickers -----
    def on_day0_select(self, _):
        dlg = wx.FileDialog(self, "Select Day 0 files or folders", style=wx.FD_OPEN | wx.FD_MULTIPLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.day0_paths = dlg.GetPaths()
            file_count = len(resolve_uploaded_paths(self.day0_paths))
            self.day0_status.SetLabel(f"✅ {file_count} Day 0 files selected.")
        else:
            self.day0_status.SetLabel("⚠️ No Day 0 files selected.")
        dlg.Destroy()

    def on_day42_select(self, _):
        dlg = wx.FileDialog(self, "Select Day 42 files or folders", style=wx.FD_OPEN | wx.FD_MULTIPLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.day42_paths = dlg.GetPaths()
            file_count = len(resolve_uploaded_paths(self.day42_paths))
            self.day42_status.SetLabel(f"✅ {file_count} Day 42 files selected.")
        else:
            self.day42_status.SetLabel("⚠️ No Day 42 files selected.")
        dlg.Destroy()

    def on_browse_csv(self, _):
        dlg = wx.FileDialog(self, "Save CSV Output As", wildcard="CSV files (*.csv)|*.csv",
                            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            self.csv_path_ctrl.SetValue(dlg.GetPath())
        dlg.Destroy()

    def on_browse_pdf(self, _):
        dlg = wx.FileDialog(self, "Save PDF Output As", wildcard="PDF files (*.pdf)|*.pdf",
                            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            self.pdf_path_ctrl.SetValue(dlg.GetPath())
        dlg.Destroy()
        
    def on_day0_select2(self, _):
        dlg = wx.FileDialog(self, "Select Day 0 files/folders (Set 2)", style=wx.FD_OPEN | wx.FD_MULTIPLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.day0_paths2 = dlg.GetPaths()
            file_count = len(resolve_uploaded_paths(self.day0_paths2))
            self.day0_status2.SetLabel(f"✅ {file_count} Day 0 files selected (Set 2).")
        else:
            self.day0_status2.SetLabel("⚠️ No Day 0 files selected (Set 2).")
        dlg.Destroy()

    def on_day42_select2(self, _):
        dlg = wx.FileDialog(self, "Select Day 42 files/folders (Set 2)", style=wx.FD_OPEN | wx.FD_MULTIPLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.day42_paths2 = dlg.GetPaths()
            file_count = len(resolve_uploaded_paths(self.day42_paths2))
            self.day42_status2.SetLabel(f"✅ {file_count} Day 42 files selected (Set 2).")
        else:
            self.day42_status2.SetLabel("⚠️ No Day 42 files selected (Set 2).")
        dlg.Destroy()

    # ----- Core actions -----
    def on_normalize(self, _):
        try:
            # --- Set 1 ---
            day0_files = resolve_uploaded_paths(self.day0_paths)
            day42_files = resolve_uploaded_paths(self.day42_paths)
            if not day0_files or not day42_files:
                self.status.SetLabel("⚠️ One or both file groups missing (Set 1).")
                return

            df_day0 = read_and_normalize_files(day0_files, "Day0")
            df_day42 = read_and_normalize_files(day42_files, "Day42")

            combined = pd.merge(df_day0, df_day42, on="repeat length", how="inner")
            day0_cols = [c for c in combined.columns if c.startswith("Day0_")]
            day42_cols = [c for c in combined.columns if c.startswith("Day42_")]

            combined["Day0_avg"] = combined[day0_cols].mean(axis=1)
            combined["Day42_avg"] = combined[day42_cols].mean(axis=1)
            combined["delta"] = combined["Day42_avg"] - combined["Day0_avg"]
            cols = [c for c in combined.columns if c != "delta"]
            combined = combined[cols + ["delta"]]

            csv_path = self.csv_path_ctrl.GetValue()
            combined.to_csv(csv_path, index=False)
            self.merged_df = combined
            msg = f"✅ Set 1 CSV saved to: {csv_path}"

            # --- Optional Set 2 ---
            day0_files2 = resolve_uploaded_paths(self.day0_paths2)
            day42_files2 = resolve_uploaded_paths(self.day42_paths2)
            self.merged_df2 = None

            if day0_files2 and day42_files2:
                df_day0_2 = read_and_normalize_files(day0_files2, "Day0")
                df_day42_2 = read_and_normalize_files(day42_files2, "Day42")

                combined2 = pd.merge(df_day0_2, df_day42_2, on="repeat length", how="inner")
                d0c2 = [c for c in combined2.columns if c.startswith("Day0_")]
                d42c2 = [c for c in combined2.columns if c.startswith("Day42_")]
                combined2["Day0_avg"] = combined2[d0c2].mean(axis=1)
                combined2["Day42_avg"] = combined2[d42c2].mean(axis=1)
                combined2["delta"] = combined2["Day42_avg"] - combined2["Day0_avg"]
                cols2 = [c for c in combined2.columns if c != "delta"]
                combined2 = combined2[cols2 + ["delta"]]
                # Save next to set1 csv with suffix
                root, ext = os.path.splitext(csv_path)
                csv2 = f"{root}_set2{ext or '.csv'}"
                combined2.to_csv(csv2, index=False)
                self.merged_df2 = combined2
                msg += f" | ✅ Set 2 CSV saved to: {csv2}"
            else:
                msg += " | (Set 2 not provided)"

            self.status.SetLabel(msg)

        except Exception as e:
            self.status.SetLabel(f"❌ Error: {e}")


    def on_plot(self, _):
        if self.merged_df is None:
            self.status.SetLabel("⚠️ No data to plot. Please normalise first.")
            return
        bin_size = self.bin_spin.GetValue()
        title = self.title_txt.GetValue()
        color = self.color_choice.GetStringSelection()
        error_mode = self.error_bar_choice.GetStringSelection()
        pdf_path = self.pdf_path_ctrl.GetValue()

        # pick a contrasting default for set 2
        color2 = "black"

        generate_delta_plot(self.merged_df, output_path=pdf_path, bin_size=bin_size,
                           title=title, color=color, error_mode=error_mode,
                           dataframe2=self.merged_df2, color2=color2)
        self.status.SetLabel(f"✅ Delta plot saved to: {pdf_path}")

    def on_view_plot(self, _event):
        if self.merged_df is None:
            self.status.SetLabel("⚠️ Please normalize data first.")
            return

        bin_size   = self.bin_spin.GetValue()
        title      = self.title_txt.GetValue()
        color      = self.color_choice.GetStringSelection()
        error_mode = self.error_bar_choice.GetStringSelection()

        # Clear panel BEFORE adding new canvas
        self.plot_panel.DestroyChildren()
        if getattr(self, "plot_canvas", None):
            self.plot_canvas.Destroy()
            self.plot_canvas = None

        # Build figure
        fig = Figure(figsize=(6, 3))
        ax  = fig.add_subplot(111)
        ax.axhline(0, color='gray', linestyle='--', linewidth=1)

        # ---- Set 1 ----
        df1   = self.merged_df.copy().sort_values("repeat length")
        df1_f = df1.iloc[::bin_size, :]
        x1, y1 = df1_f["repeat length"], df1_f["delta"]

        yerr1 = None
        if error_mode != "None":
            d0c = [c for c in df1.columns if c.startswith("Day0_")]
            d42c = [c for c in df1.columns if c.startswith("Day42_")]
            if d0c and d42c:
                d0v = df1[d0c].iloc[::bin_size, :]
                d42v = df1[d42c].iloc[::bin_size, :]
                if error_mode.startswith("Standard Deviation"):
                    yerr1 = np.sqrt(d0v.std(axis=1)**2 + d42v.std(axis=1)**2)
                else:
                    yerr1 = np.sqrt(d0v.sem(axis=1)**2 + d42v.sem(axis=1)**2)

        if yerr1 is not None:
            ax.errorbar(x1, y1, yerr=yerr1, fmt='o-', color=color, capsize=3, markersize=3, label="Set 1")
        else:
            ax.plot(x1, y1, marker='o', linestyle='-', color=color, markersize=3, label="Set 1")

        # ---- Optional Set 2 ----
        if getattr(self, "merged_df2", None) is not None:
            df2   = self.merged_df2.copy().sort_values("repeat length")
            df2_f = df2.iloc[::bin_size, :]
            x2, y2 = df2_f["repeat length"], df2_f["delta"]

            yerr2 = None
            if error_mode != "None":
                d0c2 = [c for c in df2.columns if c.startswith("Day0_")]
                d42c2 = [c for c in df2.columns if c.startswith("Day42_")]
                if d0c2 and d42c2:
                    d0v2 = df2[d0c2].iloc[::bin_size, :]
                    d42v2 = df2[d42c2].iloc[::bin_size, :]
                    if error_mode.startswith("Standard Deviation"):
                        yerr2 = np.sqrt(d0v2.std(axis=1)**2 + d42v2.std(axis=1)**2)
                    else:
                        yerr2 = np.sqrt(d0v2.sem(axis=1)**2 + d42v2.sem(axis=1)**2)

            if yerr2 is not None:
                ax.errorbar(x2, y2, yerr=yerr2, fmt='s--', color="black", capsize=3, markersize=3, label="Set 2")
            else:
                ax.plot(x2, y2, marker='s', linestyle='--', color="black", markersize=3, label="Set 2")

        ax.set_xlim(0, 200)
        ax.set_title(title)
        ax.set_xlabel("Repeat Length")
        ax.set_ylabel("Delta (Day42_avg - Day0_avg)")
        ax.legend()

        # Embed canvas
        self.plot_canvas = FigureCanvas(self.plot_panel, -1, fig)
        self.plot_canvas.draw()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.plot_canvas, 1, wx.EXPAND)
        self.plot_panel.SetSizer(sizer)
        self.plot_panel.Layout()
        self.plot_panel.Refresh()
        self.plot_panel.Update()

        wx.CallAfter(self.status.SetLabel, "📊 Delta plot displayed below.")


# -------------------------
# Main frame with Notebook
# -------------------------
class MainFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="RD and Delta plot GUI", size=(1000, 750))
        nb = wx.Notebook(self)
        self.gen_panel = GenerateHistogramsPanel(nb)
        self.analyze_panel = AnalyzeHistogramsPanel(nb)
        nb.AddPage(self.gen_panel, "Repeat Detector")
        nb.AddPage(self.analyze_panel, "Plot Delta")
        self.Centre()
        self.Show()

class App(wx.App):
    def OnInit(self):
        MainFrame()
        return True

if __name__ == "__main__":
    App(False).MainLoop()
