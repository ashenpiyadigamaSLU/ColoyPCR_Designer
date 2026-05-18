import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

class PrimerDesignerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Colony PCR Primer Designer")
        self.root.geometry("950x750")
        
        style = ttk.Style()
        style.theme_use('clam')
        
        # --- Left Panel: Inputs ---
        input_frame = ttk.LabelFrame(root, text=" Design Parameters ", padding=15)
        input_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        
        ttk.Label(input_frame, text="Existing Primers (One per line):").pack(anchor=tk.W)
        self.primer_input = scrolledtext.ScrolledText(input_frame, width=35, height=6, font=('Courier', 10))
        self.primer_input.pack(fill=tk.X, pady=5)
        self.primer_input.insert(tk.END, "GCTCGGTTTCCGACCT\nCGTTCGATCGATCGAA")
        
        ttk.Label(input_frame, text="Existing Primers Orientation:").pack(anchor=tk.W, pady=(5,0))
        self.orientation_var = tk.StringVar(value="Sense")
        ttk.Radiobutton(input_frame, text="Sense (Forces Antisense design)", variable=self.orientation_var, value="Sense").pack(anchor=tk.W)
        ttk.Radiobutton(input_frame, text="Antisense (Forces Sense design)", variable=self.orientation_var, value="Antisense").pack(anchor=tk.W)
        
        ttk.Label(input_frame, text="Target Sequence (5' -> 3'):").pack(anchor=tk.W, pady=(15,0))
        self.sequence_input = scrolledtext.ScrolledText(input_frame, width=35, height=15, font=('Courier', 10), wrap=tk.CHAR)
        self.sequence_input.pack(fill=tk.BOTH, expand=True, pady=5)
        self.sequence_input.insert(tk.END, "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGCTCGGTTTCCGACCTATTTGGCCAAATTTAAAACCCGGG")
        
        # Configure custom highlight tag for the text box
        self.sequence_input.tag_configure("highlight", background="yellow", foreground="black")
        
        self.calc_btn = ttk.Button(input_frame, text="Design Compatible Primers", command=self.design_primers)
        self.calc_btn.pack(fill=tk.X, pady=10)

        # --- Right Panel: Results ---
        result_frame = ttk.LabelFrame(root, text=" Compatible Primer Options ", padding=15)
        result_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # NEB Verification Warning Header (Fixed)
        warning_banner = ttk.Label(
            result_frame, 
            text="⚠️ WARNING: Biopython default values used.\nPlease verify the final annealing temperature on the official NEB Tm Calculator website!",
            font=('Helvetica', 9, 'bold'), padding=8, justify=tk.LEFT, anchor=tk.W
        )
        style.configure("Warning.TLabel", foreground="darkred", background="#ffe6e6")
        warning_banner.configure(style="Warning.TLabel")
        warning_banner.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(result_frame, text="💡 Click a row to highlight it in the sequence | Double-click to copy.", foreground="gray").pack(anchor=tk.W, pady=(0, 5))
        
        columns = ('seq', 'len', 'gc', 'tm', 'ta_range')
        self.tree = ttk.Treeview(result_frame, columns=columns, show='headings')
        self.tree.heading('seq', text='Sequence (5\'->3\')')
        self.tree.heading('len', text='Len')
        self.tree.heading('gc', text='GC%')
        self.tree.heading('tm', text='Tm (°C)')
        self.tree.heading('ta_range', text='Est. Ta (°C)')
        
        self.tree.column('seq', width=240, anchor=tk.W)
        self.tree.column('len', width=45, anchor=tk.CENTER)
        self.tree.column('gc', width=55, anchor=tk.CENTER)
        self.tree.column('tm', width=65, anchor=tk.CENTER)
        self.tree.column('ta_range', width=100, anchor=tk.CENTER)
        
        scrollbar = ttk.Scrollbar(result_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscroll=scrollbar.set)
        
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Table Bindings
        self.tree.bind("<<TreeviewSelect>>", self.highlight_sequence)
        self.tree.bind("<Double-1>", self.copy_sequence)
        
        # --- Bottom Status Bar ---
        self.status_bar = ttk.Label(root, text="Ready", relief=tk.SUNKEN, anchor=tk.W, padding=(5, 2))
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    def copy_sequence(self, event):
        selected_item = self.tree.selection()
        if selected_item:
            row_values = self.tree.item(selected_item[0], 'values')
            if row_values:
                primer_seq = row_values[0]
                self.root.clipboard_clear()
                self.root.clipboard_append(primer_seq)
                self.root.update()
                self.status_bar.config(text=f"📋 Copied to clipboard: {primer_seq}")

    def calculate_biopython_tm(self, seq_str):
        """Standard Biopython baseline Nearest-Neighbor calculation"""
        return mt.Tm_NN(
            Seq(seq_str), 
            dnac1=1000, 
            dnac2=0, 
            Na=50, 
            Mg=1.5, 
            dNTPs=0.2, 
            nn_table=mt.DNA_NN3
        )

    def get_gc_content(self, seq):
        g_c = seq.upper().count('G') + seq.upper().count('C')
        return (g_c / len(seq)) * 100

    def highlight_sequence(self, event):
        """Finds and highlights the selected primer (or its binding site) in the sequence box."""
        self.sequence_input.tag_remove("highlight", "1.0", tk.END)
        
        selected_item = self.tree.selection()
        if not selected_item:
            return
            
        row_values = self.tree.item(selected_item[0], 'values')
        if not row_values:
            return
            
        primer_seq = row_values[0]
        full_text = self.sequence_input.get("1.0", tk.END).upper()
        
        start_idx = full_text.find(primer_seq)
        
        if start_idx == -1:
            rev_comp_seq = str(Seq(primer_seq).reverse_complement())
            start_idx = full_text.find(rev_comp_seq)
            
        if start_idx != -1:
            lines_before = full_text[:start_idx].count('\n')
            
            if lines_before == 0:
                char_idx = start_idx
            else:
                last_newline = full_text[:start_idx].rfind('\n')
                char_idx = start_idx - last_newline - 1
                
            length = len(primer_seq)
            
            tk_start = f"{lines_before + 1}.{char_idx}"
            tk_end = f"{lines_before + 1}.{char_idx + length}"
            
            self.sequence_input.tag_add("highlight", tk_start, tk_end)
            self.sequence_input.see(tk_start)
            self.status_bar.config(text=f"Highlighted sequence position: {tk_start} to {tk_end}")
        else:
            self.status_bar.config(text="Could not map exact location on target box.")

    def design_primers(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.sequence_input.tag_remove("highlight", "1.0", tk.END)
        self.status_bar.config(text="Processing...")
            
        raw_primers = self.primer_input.get("1.0", tk.END).strip().split('\n')
        existing_primers = [p.strip().upper() for p in raw_primers if p.strip()]
        
        if not existing_primers:
            messagebox.showerror("Error", "Please enter at least one existing primer.")
            return
            
        try:
            existing_tms = [self.calculate_biopython_tm(p) for p in existing_primers]
        except Exception as e:
            messagebox.showerror("Calculation Error", f"Invalid character in existing primers.\n{e}")
            return
            
        min_exist_tm = min(existing_tms)
        max_exist_tm = max(existing_tms)
        
        if max_exist_tm - min_exist_tm > 4.0:
            messagebox.showwarning("Warning", "Your existing primers have a Tm spread > 4°C.")

        target_seq = "".join(self.sequence_input.get("1.0", tk.END).split()).upper()
        if len(target_seq) < 24:
            messagebox.showerror("Error", "Target sequence is too short.")
            return

        orientation = self.orientation_var.get()
        working_seq = target_seq if orientation == "Antisense" else str(Seq(target_seq).reverse_complement())

        valid_candidates = []
        
        for length in range(18, 25):
            for i in range(len(working_seq) - length + 1):
                candidate = working_seq[i:i+length]
                
                gc = self.get_gc_content(candidate)
                if not (40 <= gc <= 60):
                    continue
                    
                three_prime_end = candidate[-5:]
                if (three_prime_end.count('G') + three_prime_end.count('C')) > 3:
                    continue
                
                tm = self.calculate_biopython_tm(candidate)
                
                if abs(tm - min_exist_tm) <= 2.0 and abs(tm - max_exist_tm) <= 2.0:
                    pair_min_tm = min(tm, min_exist_tm)
                    ta = round(pair_min_tm - 5, 1)
                    valid_candidates.append((candidate, length, round(gc, 1), round(tm, 1), f"{ta} °C"))

        target_ideal_tm = (min_exist_tm + max_exist_tm) / 2
        valid_candidates.sort(key=lambda x: abs(x[3] - target_ideal_tm))

        if not valid_candidates:
            self.status_bar.config(text="No valid primers found.")
            messagebox.showinfo("No Results", "No matching primers found.")
        else:
            for cand in valid_candidates:
                self.tree.insert('', tk.END, values=cand)
            self.status_bar.config(text=f"Found {len(valid_candidates)} options. Click any row to map & highlight sequence.")

if __name__ == "__main__":
    root = tk.Tk()
    app = PrimerDesignerApp(root)
    root.mainloop()