# Colony PCR Primer Designer 🧬

A lightweight, standalone desktop application built in Python using **Tkinter** and **Biopython**. This tool simplifies colony PCR assay design when you already have a set of verified or universal primers in your lab and only need to design the missing matching pair (Sense or Antisense).

The core algorithm scans your template sequence and filters out candidates that do not thermodynamically align with **all** existing primers in your database, saving you bench time and master mix.

---

## 🚀 Key Features

* **Universal Compatibility Check:** Ensures designed primers match within a tight $\pm 2^\circ\text{C}$ $T_m$ window of *every* existing primer in your entered list.
* **Smart Orientation Handling:** Input your existing database as either Sense or Antisense. The app automatically handles the reverse-complementation logic to find the correct binding companion.
* **Interactive Sequence Mapping:** Clicking any designed primer option in the results table instantly highlights its exact physical binding location on your input target sequence box.
* **Colony PCR Filtering Rules:** Rejects candidates with bad 3' GC clamps (>3 G/C bases in the last 5 nucleotides) to reduce non-specific binding during rough colony lysates.
* **Quick Copy:** Double-click any row in the results table to instantly copy the primer sequence to your system clipboard.

---

## 🛠️ Installation

### Prerequisites
Make sure you have Python 3.8+ installed. You will need the `biopython` library for thermodynamic calculations.

```bash
pip install biopython
