import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

class PrimerDesignerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Colony PCR Primer Suite")
        self.root.geometry("1100x850")
        
        style = ttk.Style()
        style.theme_use('clam')
        style.configure("Warning.TLabel", foreground="darkred", background="#ffe6e6")
        
        # --- Create Master Notebook Tabs ---
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.tab1 = ttk.Frame(self.notebook)
        self.tab2 = ttk.Frame(self.notebook)
        
        self.notebook.add(self.tab1, text=" Consensus Designer (Match Database) ")
        self.notebook.add(self.tab2, text=" De Novo Pair Designer (Two Sequences) ")
        
        # --- Build Content Spaces ---
        self.setup_tab1()
        self.setup_tab2()
        
        # --- Global Bottom Status Bar ---
        self.status_bar = ttk.Label(root, text="System Ready", relief=tk.SUNKEN, anchor=tk.W, padding=(5, 2))
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    def calculate_biopython_tm(self, seq_str):
        """Standard Biopython baseline Nearest-Neighbor calculation"""
        return mt.Tm_NN(Seq(seq_str), dnac1=1000, dnac2=0, Na=50, Mg=1.5, dNTPs=0.2, nn_table=mt.DNA_NN3)

    def get_gc_content(self, seq):
        g_c = seq.upper().count('G') + seq.upper().count('C')
        return (g_c / len(seq)) * 100

    def copy_sequence_generic(self, tree):
        selected_item = tree.selection()
        if selected_item:
            row_values = tree.item(selected_item[0], 'values')
            if row_values:
                primer_seq = row_values[0] # Grab first column item
                self.root.clipboard_clear()
                self.root.clipboard_append(primer_seq)
                self.root.update()
                self.status_bar.config(text=f"📋 Copied to clipboard: {primer_seq}")

    # =========================================================================
    # TAB 1: ORIGINAL CONSENSUS SYSTEM
    # =========================================================================
    def setup_tab1(self):
        input_frame = ttk.LabelFrame(self.tab1, text=" Design Parameters & Thresholds ", padding=15)
        input_frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10, pady=10)
        
        ttk.Label(input_frame, text="Existing Primers (One per line):").pack(anchor=tk.W)
        self.primer_input = scrolledtext.ScrolledText(input_frame, width=35, height=5, font=('Courier', 10))
        self.primer_input.pack(fill=tk.X, pady=5)
        self.primer_input.insert(tk.END, "GCTCGGTTTCCGACCT\nCGTTCGATCGATCGAA")
        
        ttk.Label(input_frame, text="Existing Primers Orientation:").pack(anchor=tk.W, pady=(5,0))
        self.orientation_var = tk.StringVar(value="Sense")
        ttk.Radiobutton(input_frame, text="Sense (Forces Antisense design)", variable=self.orientation_var, value="Sense").pack(anchor=tk.W)
        ttk.Radiobutton(input_frame, text="Antisense (Forces Sense design)", variable=self.orientation_var, value="Antisense").pack(anchor=tk.W)
        
        ttk.Separator(input_frame, orient='horizontal').pack(fill='x', pady=10)
        ttk.Label(input_frame, text="🔧 Settings Engine", font=('Helvetica', 10, 'bold')).pack(anchor=tk.W, pady=(0,5))
        
        len_frame = ttk.Frame(input_frame)
        len_frame.pack(fill=tk.X, pady=2)
        ttk.Label(len_frame, text="Length Range (bp):").pack(side=tk.LEFT)
        self.min_len_entry = ttk.Entry(len_frame, width=5)
        self.min_len_entry.pack(side=tk.LEFT, padx=2)
        self.min_len_entry.insert(0, "18")
        ttk.Label(len_frame, text="to").pack(side=tk.LEFT)
        self.max_len_entry = ttk.Entry(len_frame, width=5)
        self.max_len_entry.pack(side=tk.LEFT, padx=2)
        self.max_len_entry.insert(0, "24")
        
        ttk.Label(input_frame, text="Min GC Content (%):").pack(anchor=tk.W, pady=(5,0))
        self.min_gc_val = tk.IntVar(value=40)
        self.min_gc_lbl = ttk.Label(input_frame, text="40%")
        self.min_gc_lbl.pack(anchor=tk.E)
        ttk.Scale(input_frame, from_=20, to=80, variable=self.min_gc_val, command=lambda e: self.min_gc_lbl.config(text=f"{self.min_gc_val.get()}%")).pack(fill=tk.X)
        
        ttk.Label(input_frame, text="Max GC Content (%):").pack(anchor=tk.W, pady=(5,0))
        self.max_gc_val = tk.IntVar(value=60)
        self.max_gc_lbl = ttk.Label(input_frame, text="60%")
        self.max_gc_lbl.pack(anchor=tk.E)
        ttk.Scale(input_frame, from_=40, to=100, variable=self.max_gc_val, command=lambda e: self.max_gc_lbl.config(text=f"{self.max_gc_val.get()}%")).pack(fill=tk.X)
        
        ctrl_grid = ttk.Frame(input_frame)
        ctrl_grid.pack(fill=tk.X, pady=10)
        ttk.Label(ctrl_grid, text="Max G/C in last 5bp (Clamp):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.max_clamp_entry = ttk.Entry(ctrl_grid, width=5)
        self.max_clamp_entry.grid(row=0, column=1, sticky=tk.E, pady=2)
        self.max_clamp_entry.insert(0, "3")
        
        ttk.Label(ctrl_grid, text="Max Tm Variance Allowed (°C):").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.max_var_entry = ttk.Entry(ctrl_grid, width=5)
        self.max_var_entry.grid(row=1, column=1, sticky=tk.E, pady=2)
        self.max_var_entry.insert(0, "2.0")
        
        ttk.Label(input_frame, text="Target Sequence (5' -> 3'):").pack(anchor=tk.W, pady=(5,0))
        self.sequence_input = scrolledtext.ScrolledText(input_frame, width=35, height=8, font=('Courier', 10), wrap=tk.CHAR)
        self.sequence_input.pack(fill=tk.BOTH, expand=True, pady=5)
        self.sequence_input.insert(tk.END, "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGCTCGGTTTCCGACCTATTTGGCCAAATTTAAAACCCGGG")
        self.sequence_input.tag_configure("highlight", background="yellow", foreground="black")
        
        ttk.Button(input_frame, text="Design Compatible Primers", command=self.design_primers_tab1).pack(fill=tk.X, pady=10)

        # Result Panel Right Side
        result_frame = ttk.LabelFrame(self.tab1, text=" Compatible Primer Options ", padding=15)
        result_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        warning_banner = ttk.Label(result_frame, text="⚠️ WARNING: Biopython default values used. Verify temperatures on NEB website!", font=('Helvetica', 9, 'bold'), padding=8, justify=tk.LEFT, anchor=tk.W)
        warning_banner.configure(style="Warning.TLabel")
        warning_banner.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(result_frame, text="💡 Click a row to highlight it in the sequence | Double-click to copy.", foreground="gray").pack(anchor=tk.W, pady=(0, 5))
        
        columns = ('seq', 'len', 'gc', 'tm', 'ta_range')
        self.tree1 = ttk.Treeview(result_frame, columns=columns, show='headings')
        for col in columns: self.tree1.heading(col, text=col.upper())
        self.tree1.column('seq', width=240, anchor=tk.W)
        self.tree1.column('len', width=45, anchor=tk.CENTER)
        self.tree1.column('gc', width=55, anchor=tk.CENTER)
        self.tree1.column('tm', width=65, anchor=tk.CENTER)
        self.tree1.column('ta_range', width=100, anchor=tk.CENTER)
        
        scrollbar = ttk.Scrollbar(result_frame, orient=tk.VERTICAL, command=self.tree1.yview)
        self.tree1.configure(yscroll=scrollbar.set)
        self.tree1.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.tree1.bind("<<TreeviewSelect>>", self.highlight_sequence_tab1)
        self.tree1.bind("<Double-1>", lambda e: self.copy_sequence_generic(self.tree1))

    def highlight_sequence_tab1(self, event):
        self.sequence_input.tag_remove("highlight", "1.0", tk.END)
        selected = self.tree1.selection()
        if not selected: return
        row_values = self.tree1.item(selected[0], 'values')
        if not row_values: return
        
        primer_seq = row_values[0]
        full_text = self.sequence_input.get("1.0", tk.END).upper()
        
        start_idx = full_text.find(primer_seq)
        if start_idx == -1:
            start_idx = full_text.find(str(Seq(primer_seq).reverse_complement()))
            
        if start_idx != -1:
            lines = full_text[:start_idx].count('\n')
            char_idx = start_idx if lines == 0 else start_idx - full_text[:start_idx].rfind('\n') - 1
            tk_start = f"{lines + 1}.{char_idx}"
            tk_end = f"{lines + 1}.{char_idx + len(primer_seq)}"
            self.sequence_input.tag_add("highlight", tk_start, tk_end)
            self.sequence_input.see(tk_start)

    def design_primers_tab1(self):
        for item in self.tree1.get_children(): self.tree1.delete(item)
        try:
            min_len, max_len = int(self.min_len_entry.get()), int(self.max_len_entry.get())
            min_gc, max_gc = self.min_gc_val.get(), self.max_gc_val.get()
            max_clamp, max_var = int(self.max_clamp_entry.get()), float(self.max_var_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Check your parameter criteria entries.")
            return
            
        raw_primers = self.primer_input.get("1.0", tk.END).strip().split('\n')
        existing_primers = [p.strip().upper() for p in raw_primers if p.strip()]
        if not existing_primers: return
        
        existing_tms = [self.calculate_biopython_tm(p) for p in existing_primers]
        min_exist_tm, max_exist_tm = min(existing_tms), max(existing_tms)
        target_seq = "".join(self.sequence_input.get("1.0", tk.END).split()).upper()
        
        orientation = self.orientation_var.get()
        working_seq = target_seq if orientation == "Antisense" else str(Seq(target_seq).reverse_complement())
        valid_candidates = []
        
        for length in range(min_len, max_len + 1):
            for i in range(len(working_seq) - length + 1):
                candidate = working_seq[i:i+length]
                gc = self.get_gc_content(candidate)
                if not (min_gc <= gc <= max_gc): continue
                if (candidate[-5:].count('G') + candidate[-5:].count('C')) > max_clamp: continue
                tm = self.calculate_biopython_tm(candidate)
                
                if abs(tm - min_exist_tm) <= max_var and abs(tm - max_exist_tm) <= max_var:
                    ta = round(min(tm, min_exist_tm) - 5, 1)
                    valid_candidates.append((candidate, length, round(gc, 1), round(tm, 1), f"{ta} °C"))
                    
        valid_candidates.sort(key=lambda x: abs(x[3] - (min_exist_tm + max_exist_tm)/2))
        for cand in valid_candidates: self.tree1.insert('', tk.END, values=cand)


    # =========================================================================
    # TAB 2: DE NOVO PRIMER PAIR GENERATOR
    # =========================================================================
    def setup_tab2(self):
        # Left Input Control Frame
        input_frame = ttk.LabelFrame(self.tab2, text=" Sequence & Dynamic Pair Parameters ", padding=15)
        input_frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10, pady=10)
        
        # Two input areas for sequence frames
        ttk.Label(input_frame, text="Forward Target Template Sequence (5' -> 3'):").pack(anchor=tk.W)
        self.fwd_seq_input = scrolledtext.ScrolledText(input_frame, width=38, height=8, font=('Courier', 10), wrap=tk.CHAR)
        self.fwd_seq_input.pack(fill=tk.X, pady=5)
        self.fwd_seq_input.insert(tk.END, "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC")
        self.fwd_seq_input.tag_configure("highlight", background="lightgreen", foreground="black")
        
        ttk.Label(input_frame, text="Reverse Target Template Sequence (5' -> 3'):").pack(anchor=tk.W, pady=(5,0))
        self.rev_seq_input = scrolledtext.ScrolledText(input_frame, width=38, height=8, font=('Courier', 10), wrap=tk.CHAR)
        self.rev_seq_input.pack(fill=tk.X, pady=5)
        self.rev_seq_input.insert(tk.END, "GCTCGGTTTCCGACCTATTTGGCCAAATTTAAAACCCGGG")
        self.rev_seq_input.tag_configure("highlight", background="lightpink", foreground="black")
        
        ttk.Separator(input_frame, orient='horizontal').pack(fill='x', pady=10)
        
        # Dual-Pair Threshold Controls
        ttk.Label(input_frame, text="⚙️ Pair Criteria Limits", font=('Helvetica', 10, 'bold')).pack(anchor=tk.W, pady=(0,5))
        
        len_frame = ttk.Frame(input_frame)
        len_frame.pack(fill=tk.X, pady=2)
        ttk.Label(len_frame, text="Length Bounds (bp):").pack(side=tk.LEFT)
        self.tab2_min_len = ttk.Entry(len_frame, width=5)
        self.tab2_min_len.pack(side=tk.LEFT, padx=2)
        self.tab2_min_len.insert(0, "18")
        ttk.Label(len_frame, text="to").pack(side=tk.LEFT)
        self.tab2_max_len = ttk.Entry(len_frame, width=5)
        self.tab2_max_len.pack(side=tk.LEFT, padx=2)
        self.tab2_max_len.insert(0, "24")
        
        ttk.Label(input_frame, text="Acceptable GC Content window (%):").pack(anchor=tk.W, pady=(5,0))
        gc_bounds_frame = ttk.Frame(input_frame)
        gc_bounds_frame.pack(fill=tk.X)
        self.tab2_min_gc = ttk.Entry(gc_bounds_frame, width=5)
        self.tab2_min_gc.pack(side=tk.LEFT)
        self.tab2_min_gc.insert(0, "40")
        ttk.Label(gc_bounds_frame, text="%  to").pack(side=tk.LEFT, padx=4)
        self.tab2_max_gc = ttk.Entry(gc_bounds_frame, width=5)
        self.tab2_max_gc.pack(side=tk.LEFT)
        self.tab2_max_gc.insert(0, "60")
        ttk.Label(gc_bounds_frame, text="%").pack(side=tk.LEFT)
        
        grid_frame = ttk.Frame(input_frame)
        grid_frame.pack(fill=tk.X, pady=10)
        ttk.Label(grid_frame, text="Max 3' End G/C Clamp:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.tab2_clamp = ttk.Entry(grid_frame, width=5)
        self.tab2_clamp.grid(row=0, column=1, sticky=tk.E, pady=2)
        self.tab2_clamp.insert(0, "3")
        
        ttk.Label(grid_frame, text="Max Delta-Tm Between Pairs (°C):").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.tab2_delta_tm = ttk.Entry(grid_frame, width=5)
        self.tab2_delta_tm.grid(row=1, column=1, sticky=tk.E, pady=2)
        self.tab2_delta_tm.insert(0, "1.5")
        
        ttk.Button(input_frame, text="Design Matching Primer Pairs", command=self.design_pairs_tab2).pack(fill=tk.X, pady=10)

        # Right Result Panel Space
        result_frame = ttk.LabelFrame(self.tab2, text=" Discovered Primer Pairs Match Lists ", padding=15)
        result_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        warning_banner = ttk.Label(result_frame, text="⚠️ WARNING: Biopython default values used. Verify temperatures on NEB website!", font=('Helvetica', 9, 'bold'), padding=8, justify=tk.LEFT, anchor=tk.W)
        warning_banner.configure(style="Warning.TLabel")
        warning_banner.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(result_frame, text="💡 Click a row to map BOTH primers | Double-click row to copy Forward primer.", foreground="gray").pack(anchor=tk.W, pady=(0, 5))
        
        columns = ('fwd_seq', 'fwd_tm', 'rev_seq', 'rev_tm', 'delta_tm', 'ta_value')
        self.tree2 = ttk.Treeview(result_frame, columns=columns, show='headings')
        
        self.tree2.heading('fwd_seq', text='Forward Primer')
        self.tree2.heading('fwd_tm', text='F-Tm (°C)')
        self.tree2.heading('rev_seq', text='Reverse Primer')
        self.tree2.heading('rev_tm', text='R-Tm (°C)')
        self.tree2.heading('delta_tm', text='ΔTm')
        self.tree2.heading('ta_value', text='Pair Ta')
        
        self.tree2.column('fwd_seq', width=160, anchor=tk.W)
        self.tree2.column('fwd_tm', width=55, anchor=tk.CENTER)
        self.tree2.column('rev_seq', width=160, anchor=tk.W)
        self.tree2.column('rev_tm', width=55, anchor=tk.CENTER)
        self.tree2.column('delta_tm', width=45, anchor=tk.CENTER)
        self.tree2.column('ta_value', width=70, anchor=tk.CENTER)
        
        scrollbar = ttk.Scrollbar(result_frame, orient=tk.VERTICAL, command=self.tree2.yview)
        self.tree2.configure(yscroll=scrollbar.set)
        self.tree2.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.tree2.bind("<<TreeviewSelect>>", self.highlight_sequence_tab2)
        self.tree2.bind("<Double-1>", lambda e: self.copy_sequence_generic(self.tree2))

    def highlight_sequence_tab2(self, event):
        # Clear highlights inside both target windows
        self.fwd_seq_input.tag_remove("highlight", "1.0", tk.END)
        self.rev_seq_input.tag_remove("highlight", "1.0", tk.END)
        
        selected = self.tree2.selection()
        if not selected: return
        row_values = self.tree2.item(selected[0], 'values')
        if not row_values: return
        
        fwd_primer, _, rev_primer, _, _, _ = row_values
        
        # 1. Map Forward Primer (Acts directly on Forward Input box strand)
        fwd_text = self.fwd_seq_input.get("1.0", tk.END).upper()
        f_start = fwd_text.find(fwd_primer)
        if f_start != -1:
            f_lines = fwd_text[:f_start].count('\n')
            f_char = f_start if f_lines == 0 else f_start - fwd_text[:f_start].rfind('\n') - 1
            self.fwd_seq_input.tag_add("highlight", f"{f_lines+1}.{f_char}", f"{f_lines+1}.{f_char+len(fwd_primer)}")
            
        # 2. Map Reverse Primer (The designed candidate is written 5'->3', so look for its complementary sequence in Reverse Box)
        rev_text = self.rev_seq_input.get("1.0", tk.END).upper()
        r_start = rev_text.find(str(Seq(rev_primer).reverse_complement()))
        if r_start != -1:
            r_lines = rev_text[:r_start].count('\n')
            r_char = r_start if r_lines == 0 else r_start - rev_text[:r_start].rfind('\n') - 1
            self.rev_seq_input.tag_add("highlight", f"{r_lines+1}.{r_char}", f"{r_lines+1}.{r_char+len(rev_primer)}")

    def design_pairs_tab2(self):
        for item in self.tree2.get_children(): self.tree2.delete(item)
        self.fwd_seq_input.tag_remove("highlight", "1.0", tk.END)
        self.rev_seq_input.tag_remove("highlight", "1.0", tk.END)
        
        try:
            min_len = int(self.tab2_min_len.get())
            max_len = int(self.tab2_max_len.get())
            min_gc = int(self.tab2_min_gc.get())
            max_gc = int(self.tab2_max_gc.get())
            max_clamp = int(self.tab2_clamp.get())
            max_delta = float(self.tab2_delta_tm.get())
        except ValueError:
            messagebox.showerror("Configuration Error", "Please verify setup input numbers.")
            return
            
        fwd_template = "".join(self.fwd_seq_input.get("1.0", tk.END).split()).upper()
        rev_template = "".join(self.rev_seq_input.get("1.0", tk.END).split()).upper()
        
        if len(fwd_template) < min_len or len(rev_template) < min_len:
            messagebox.showerror("Data Error", "Sequences entered are too short for target window lengths.")
            return
            
        # Step A: Harvest all candidate Forward Primers (from the forward strand)
        fwd_candidates = []
        for length in range(min_len, max_len + 1):
            for i in range(len(fwd_template) - length + 1):
                sub_seq = fwd_template[i:i+length]
                gc = self.get_gc_content(sub_seq)
                if not (min_gc <= gc <= max_gc): continue
                if (sub_seq[-5:].count('G') + sub_seq[-5:].count('C')) > max_clamp: continue
                tm = self.calculate_biopython_tm(sub_seq)
                fwd_candidates.append((sub_seq, tm))
                
        # Step B: Harvest all candidate Reverse Primers (take reverse-complement of the reverse template strand)
        rev_strand_working = str(Seq(rev_template).reverse_complement())
        rev_candidates = []
        for length in range(min_len, max_len + 1):
            for i in range(len(rev_strand_working) - length + 1):
                sub_seq = rev_strand_working[i:i+length]
                gc = self.get_gc_content(sub_seq)
                if not (min_gc <= gc <= max_gc): continue
                if (sub_seq[-5:].count('G') + sub_seq[-5:].count('C')) > max_clamp: continue
                tm = self.calculate_biopython_tm(sub_seq)
                rev_candidates.append((sub_seq, tm))
                
        # Step C: Pairwise Cross-Matching Analysis Matrix
        matched_pairs = []
        for f_seq, f_tm in fwd_candidates:
            for r_seq, r_tm in rev_candidates:
                diff = abs(f_tm - r_tm)
                if diff <= max_delta:
                    ta = round(min(f_tm, r_tm) - 5, 1)
                    matched_pairs.append((f_seq, round(f_tm, 1), r_seq, round(r_tm, 1), round(diff, 1), f"{ta} °C"))
                    
        # Sort by lowest Delta Tm difference between pairs
        matched_pairs.sort(key=lambda x: x[4])
        
        if not matched_pairs:
            self.status_bar.config(text="No matching compatible primer pairs found.")
            messagebox.showinfo("Finished", "Zero sequence combinations survived pair filtration parameters.")
        else:
            for pair in matched_pairs:
                self.tree2.insert('', tk.END, values=pair)
            self.status_bar.config(text=f"Found {len(matched_pairs)} optimized primer pair configurations.")

if __name__ == "__main__":
    root = tk.Tk()
    app = PrimerDesignerApp(root)
    root.mainloop()